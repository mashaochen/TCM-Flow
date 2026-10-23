# -*- coding: utf-8 -*-
"""
MD 轨迹分析 agent（MDTraj：RMSD / 回转半径 Rg / RMSF）

优先使用 MDTraj；不可用时降级为纯 numpy 实现
（手动解析 DCD + Kabsch 算法 + 标准公式），保证七阶段全流程
在没有 mdtraj 的 Windows 部署环境下也能跑通。

输出：时间序列 CSV + 汇总指标 dict
"""
import os
import struct
import numpy as np


# ---------------------------------------------------------------------------
# DCD 二进制格式读取器（标准 CHARMM/NAMD DCD 格式）
# ---------------------------------------------------------------------------
def _read_dcd(dcd_path):
    """读取 DCD 文件，返回 (n_atoms, n_frames, frames)。
    frames 是 list of np.array(n_atoms, 3)，单位 Å（DCD原始单位）。
    """
    with open(dcd_path, "rb") as f:
        # --- 第一块：头部 ---
        block_sz = struct.unpack("i", f.read(4))[0]
        magic = f.read(4)
        if magic not in (b"CORD", b"VELD"):
            raise RuntimeError(f"不是有效的 DCD 文件: magic={magic}")

        nfi = struct.unpack("i", f.read(4))[0]       # 初始步数
        nstep_save = struct.unpack("i", f.read(4))[0] # 保存间隔（步数）
        nfiav = struct.unpack("i", f.read(4))[0]      # 最后一步
        f.read(4 * 17)  # 跳过其余17个int字段
        f.read(4)  # 块尾4字节

        # --- 第二块：标题 ---
        title_sz = struct.unpack("i", f.read(4))[0]
        ntitle = struct.unpack("i", f.read(4))[0]
        f.read(title_sz - 4)  # 跳过标题内容
        f.read(4)  # 块尾

        # --- 第三块：natom ---
        natom_sz = struct.unpack("i", f.read(4))[0]
        natom = struct.unpack("i", f.read(4))[0]
        f.read(4)  # 块尾

        n_frames = (nfiav - nfi) // nstep_save + 1 if nstep_save > 0 else 0

        # --- 读取帧数据 ---
        frames = []
        coords_sz = natom * 4  # float32
        try:
            while True:
                # 检查第一块：可能是晶胞信息(48字节=6个double)，也可能直接是X坐标
                sz = struct.unpack("i", f.read(4))[0]
                if sz == 48:
                    # 晶胞块（6个double：a,b,c,alpha,beta,gamma），跳过
                    f.read(48)
                    f.read(4)  # 块尾
                    # 再读下一块，应该就是X坐标了
                    sz = struct.unpack("i", f.read(4))[0]

                if sz != coords_sz:
                    break
                x = np.frombuffer(f.read(sz), dtype=np.float32).astype(np.float64)
                f.read(4)  # 块尾
                # Y
                sz = struct.unpack("i", f.read(4))[0]
                y = np.frombuffer(f.read(sz), dtype=np.float32).astype(np.float64)
                f.read(4)
                # Z
                sz = struct.unpack("i", f.read(4))[0]
                z = np.frombuffer(f.read(sz), dtype=np.float32).astype(np.float64)
                f.read(4)

                coords = np.column_stack([x, y, z])
                frames.append(coords)
        except Exception:
            pass

        return natom, len(frames), frames


# ---------------------------------------------------------------------------
# Kabsch RMSD + Rg + RMSF（纯 numpy）
# ---------------------------------------------------------------------------
def _kabsch_rmsd(P, Q):
    """Kabsch算法，P和Q形状(n,3)，返回RMSD。"""
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign = np.ones(3)
    sign[-1] = np.sign(d)
    R = Vt.T @ np.diag(sign) @ U.T
    P_aligned = P_c @ R.T
    diff = P_aligned - Q_c
    return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))


def _compute_rg(coords):
    """回转半径（假设所有原子质量相同）。"""
    com = coords.mean(axis=0)
    diff = coords - com
    return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))


def _get_protein_atom_indices(pdb_path):
    """从PDB中获取蛋白原子（ATOM记录）的索引，对应DCD中的原子顺序。"""
    indices = []
    with open(pdb_path, encoding="utf-8") as f:
        idx = 0
        for line in f:
            if line.startswith("ATOM"):
                indices.append(idx)
                idx += 1
            elif line.startswith("HETATM"):
                idx += 1
            # 其他行忽略（TER, END等不增加原子数）
    return indices


# ---------------------------------------------------------------------------
# 降级分析入口
# ---------------------------------------------------------------------------
def _analysis_numpy(dcd_path, top_pdb, workdir, report_interval=None, dt_fs=None):
    """纯 numpy 实现的轨迹分析（mdtraj 不可用时的降级方案）。"""
    os.makedirs(workdir, exist_ok=True)

    # 读取DCD
    n_atoms_total, n_frames, frames = _read_dcd(dcd_path)
    if n_frames == 0:
        raise RuntimeError(f"无法读取 DCD 轨迹: {dcd_path}")

    # 获取蛋白原子索引
    protein_idx = _get_protein_atom_indices(top_pdb)
    n_protein = len(protein_idx)

    # 只提取蛋白原子坐标（Å -> nm）
    protein_frames = [f[protein_idx] / 10.0 for f in frames]
    ref_frame = protein_frames[0]

    # RMSD
    rmsd = np.array([_kabsch_rmsd(protein_frames[i], ref_frame)
                     for i in range(n_frames)])

    # Rg
    rg = np.array([_compute_rg(protein_frames[i]) for i in range(n_frames)])

    # RMSF（先对齐所有帧到参考结构）
    aligned_frames = []
    for i in range(n_frames):
        P = protein_frames[i] - protein_frames[i].mean(axis=0)
        Q = ref_frame - ref_frame.mean(axis=0)
        H = P.T @ Q
        U, S, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        sign = np.ones(3)
        sign[-1] = np.sign(d)
        R = Vt.T @ np.diag(sign) @ U.T
        aligned_frames.append(P @ R.T)
    aligned_stack = np.stack(aligned_frames)
    mean_pos = aligned_stack.mean(axis=0)
    diff = aligned_stack - mean_pos
    rmsf = np.sqrt(np.mean(np.sum(diff ** 2, axis=2), axis=0))

    # 时间轴（ps）
    frame_ps = 1.0
    if report_interval and dt_fs:
        frame_ps = report_interval * dt_fs / 1000.0
    time_ps = np.arange(n_frames, dtype=float) * frame_ps

    # 保存 CSV
    csv_path = os.path.join(workdir, "analysis_timeseries.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("frame,time_ps,rmsd_nm,rg_nm,rmsf_nm\n")
        for i in range(n_frames):
            f.write(f"{i},{time_ps[i]:.2f},{rmsd[i]:.4f},{rg[i]:.4f},{rmsf[i]:.4f}\n")

    metrics = {
        "n_frames": int(n_frames),
        "n_protein_atoms": int(n_protein),
        "rmsd_mean_nm": round(float(rmsd.mean()), 4),
        "rmsd_max_nm": round(float(rmsd.max()), 4),
        "rg_mean_nm": round(float(rg.mean()), 4),
        "rg_final_nm": round(float(rg[-1]), 4),
        "rmsf_mean_nm": round(float(rmsf.mean()), 4),
        "rmsf_max_nm": round(float(rmsf.max()), 4),
    }
    return {
        "metrics": metrics,
        "files": {"csv": csv_path},
        "rmsd": [round(float(x), 4) for x in rmsd],
        "rg": [round(float(x), 4) for x in rg],
        "rmsf": [round(float(x), 4) for x in rmsf],
    }


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------
def run_analysis(dcd_path, top_pdb, workdir, report_interval=None, dt_fs=None):
    """
    分析 DCD 轨迹，返回汇总指标 dict。
    report_interval / dt_fs 用于把帧序号换算成真实时间（ps）；
    缺省时退化为帧序号。
    """
    os.makedirs(workdir, exist_ok=True)

    # 优先尝试 MDTraj
    try:
        import mdtraj as md
        t = md.load(dcd_path, top=top_pdb)
        protein = t.topology.select("protein")
        if len(protein) == 0:
            protein = np.arange(t.n_atoms)

        rmsd = md.rmsd(t, t, 0, atom_indices=protein)
        rg = md.compute_rg(t.atom_slice(protein))
        rmsf = md.rmsf(t, t, 0, atom_indices=protein)

        # 时间轴（ps）
        frame_ps = 1.0
        if report_interval and dt_fs:
            frame_ps = report_interval * dt_fs / 1000.0
        time_ps = np.arange(t.n_frames, dtype=float) * frame_ps

        # 保存 CSV
        csv_path = os.path.join(workdir, "analysis_timeseries.csv")
        with open(csv_path, "w", encoding="utf-8") as f:
            f.write("frame,time_ps,rmsd_nm,rg_nm,rmsf_nm\n")
            for i in range(t.n_frames):
                f.write(f"{i},{time_ps[i]:.2f},{rmsd[i]:.4f},{rg[i]:.4f},{rmsf[i]:.4f}\n")

        metrics = {
            "n_frames": int(t.n_frames),
            "n_protein_atoms": int(len(protein)),
            "rmsd_mean_nm": round(float(rmsd.mean()), 4),
            "rmsd_max_nm": round(float(rmsd.max()), 4),
            "rg_mean_nm": round(float(rg.mean()), 4),
            "rg_final_nm": round(float(rg[-1]), 4),
            "rmsf_mean_nm": round(float(rmsf.mean()), 4),
            "rmsf_max_nm": round(float(rmsf.max()), 4),
        }
        return {
            "metrics": metrics,
            "files": {"csv": csv_path},
            "rmsd": [round(float(x), 4) for x in rmsd],
            "rg": [round(float(x), 4) for x in rg],
            "rmsf": [round(float(x), 4) for x in rmsf],
        }
    except ImportError:
        # mdtraj 不可用，降级为纯 numpy 实现
        return _analysis_numpy(dcd_path, top_pdb, workdir, report_interval, dt_fs)
