// TCM-Flow 前端逻辑
(function () {
  "use strict";

  const EXAMPLES = [
    { name: "柚皮素", smiles: "O=C1CC(c2ccc(O)cc2)Oc2cc(O)ccc12" },
    { name: "芹菜素", smiles: "O=c1cc(-c2ccc(O)cc2)oc2cc(O)cc(O)c12" },
    { name: "阿司匹林", smiles: "CC(=O)Oc1ccccc1C(=O)O" },
    { name: "白藜芦醇", smiles: "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1" },
    { name: "阿魏酸", smiles: "COc1cc(/C=C/C(=O)O)ccc1O" },
  ];
  const PDB_EXAMPLES = [
    { id: "2J0L", name: "FAK 激酶" },
    { id: "1STP", name: "链霉亲和素" },
    { id: "4LLD", name: "雌激素受体" },
  ];

  let currentUser = null;
  const $ = (s) => document.querySelector(s);
  const $$ = (s) => document.querySelectorAll(s);

  const STAGE_NAMES = { 1: "分子对接", 2: "xtb 优化", 3: "再对接", 4: "分子动力学", 5: "性质提取", 6: "轨迹分析", 7: "生成报告" };

  async function init() {
    renderChips();
    try {
      const r = await fetch("/api/auth/me"); const d = await r.json();
      if (d.logged_in) { currentUser = d; showWorkspace(); }
      else showLoginBanner();
    } catch (e) { showLoginBanner(); }
    bindEvents();
  }

  function renderChips() {
    $("#smiles-chips").innerHTML = EXAMPLES.map((e) =>
      `<span class="chip" data-smiles="${e.smiles}">${e.name} <span class="s">${e.smiles.length > 20 ? e.smiles.slice(0, 20) + "…" : e.smiles}</span></span>`).join("");
    $("#pdb-chips").innerHTML = PDB_EXAMPLES.map((e) =>
      `<span class="chip" data-pdb="${e.id}">${e.id} <span class="s">${e.name}</span></span>`).join("");
  }

  function showWorkspace() {
    $("#login-banner").style.display = "none";
    $("#workspace").style.display = "block";
    $("#nav-auth").innerHTML = `<span class="user-chip">👤 ${currentUser.username}</span>
      ${currentUser.role === "admin" ? `<a href="/admin">⚙ 管理后台</a>` : ""}
      <button onclick="logout()">退出</button>`;
  }
  function showLoginBanner() {
    $("#login-banner").style.display = "block";
    $("#workspace").style.display = "none";
    $("#nav-auth").innerHTML = `<a href="/login" class="btn-primary">登录 / 注册</a>`;
  }
  window.logout = async function () { await fetch("/api/auth/logout", { method: "POST" }); location.reload(); };

  function showNotice(id, msg, type) { const el = $(id); el.className = "notice " + (type || "error"); el.textContent = msg; }
  function hideNotice(id) { $(id).className = "notice"; }

  function bindEvents() {
    $$("#smiles-chips .chip").forEach((c) => c.addEventListener("click", () => { $("#smiles").value = c.dataset.smiles; }));
    $$("#pdb-chips .chip").forEach((c) => c.addEventListener("click", () => { $("#pdb-id").value = c.dataset.pdb; }));
    $$("#pdb-tabs button").forEach((b) => b.addEventListener("click", () => {
      $$("#pdb-tabs button").forEach((x) => x.classList.remove("active"));
      b.classList.add("active");
      const m = b.dataset.mode;
      $("#pdb-panel-pdb_id").style.display = m === "pdb_id" ? "block" : "none";
      $("#pdb-panel-upload").style.display = m === "upload" ? "block" : "none";
      $("#pdb-panel-paste").style.display = m === "paste" ? "block" : "none";
    }));
    $("#pdb-file").addEventListener("change", async (e) => {
      const f = e.target.files[0];
      if (f) { window._pdbText = await f.text(); showNotice("#pdb-notice", `已读取 ${f.name}`, "ok"); }
    });
    $("#btn-run").addEventListener("click", runPipeline);
    $("#btn-query").addEventListener("click", queryReport);
    $("#query-id").addEventListener("keydown", (e) => { if (e.key === "Enter") queryReport(); });
  }

  async function queryReport() {
    hideNotice("#query-notice");
    const id = $("#query-id").value.trim();
    if (!id) { showNotice("#query-notice", "请输入任务编号"); return; }
    try {
      const r = await fetch(`/api/full/${id}`);
      const d = await r.json();
      if (!d.success) { showNotice("#query-notice", d.error || "任务不存在"); return; }
      const job = d.job;
      if (job.status === "done") {
        renderReport(job);
        showNotice("#query-notice", `已加载任务 ${job.id} 的报告`, "ok");
      } else if (job.status === "error") {
        showNotice("#query-notice", "该任务计算失败：" + (job.error || "未知错误"));
      } else {
        const map = { queued: "排队中", running: "计算中" };
        showNotice("#query-notice", `任务 ${job.id} 状态：${map[job.status] || job.status}（尚未完成，请稍后再查）`);
      }
    } catch (e) {
      showNotice("#query-notice", "查询失败：" + e);
    }
  }

  function getProteinPayload() {
    const mode = $("#pdb-tabs button.active").dataset.mode;
    if (mode === "pdb_id") {
      const id = $("#pdb-id").value.trim().toUpperCase();
      if (!id) { showNotice("#run-notice", "请输入蛋白 PDB ID"); return null; }
      return { receptor_source: "pdb_id", pdb_id: id };
    }
    const text = (mode === "upload" ? window._pdbText : $("#pdb-text").value) || "";
    if (!text || text.length < 100) { showNotice("#run-notice", "请上传或粘贴有效 PDB"); return null; }
    return { receptor_source: "upload", pdb_text: text, pdb_id: "upload" };
  }

  async function runPipeline() {
    hideNotice("#run-notice");
    const smiles = $("#smiles").value.trim();
    if (!smiles) { showNotice("#run-notice", "请输入配体 SMILES"); return; }
    const prot = getProteinPayload();
    if (!prot) return;

    const payload = {
      smiles, ligand_name: "",
      receptor_source: prot.receptor_source, pdb_id: prot.pdb_id,
      pdb_text: prot.pdb_text || "",
      dock_params: {
        exhaustiveness: parseInt($("#p-exh").value) || 8,
        num_modes: parseInt($("#p-modes").value) || 9,
        cpu: parseInt($("#p-cpu").value) || 0,
        seed: parseInt($("#p-seed").value) || 42,
      },
    };

    $("#btn-run").disabled = true;
    $("#result-card").classList.add("result-hidden");
    $("#progress-card").classList.remove("result-hidden");
    $("#progress-msg").textContent = "任务已提交，正在排队…";

    try {
      const r = await fetch("/api/full/run", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) });
      const d = await r.json();
      if (!d.success) { showNotice("#run-notice", d.error); $("#progress-card").classList.add("result-hidden"); $("#btn-run").disabled = false; return; }
      pollJob(d.job_id);
    } catch (e) {
      showNotice("#run-notice", "提交失败：" + e);
      $("#progress-card").classList.add("result-hidden"); $("#btn-run").disabled = false;
    }
  }

  async function pollJob(jobId) {
    const timer = setInterval(async () => {
      try {
        const r = await fetch(`/api/full/${jobId}`); const d = await r.json();
        if (!d.success) { clearInterval(timer); showNotice("#run-notice", d.error); $("#btn-run").disabled = false; return; }
        const job = d.job;
        if (job.status === "queued") {
          $("#progress-msg").textContent = "任务排队中，请稍候…";
        } else if (job.status === "running") {
          const p = (job.result && job.result.progress) || {};
          const nm = STAGE_NAMES[p.stage] || "";
          $("#progress-msg").textContent = `正在执行：${nm}（${p.msg || ""}）`;
        } else if (job.status === "done") {
          clearInterval(timer);
          $("#progress-card").classList.add("result-hidden");
          renderReport(job);
          $("#btn-run").disabled = false;
        } else if (job.status === "error") {
          clearInterval(timer);
          $("#progress-card").classList.add("result-hidden");
          showNotice("#run-notice", "计算失败：" + (job.error || "未知错误"));
          $("#btn-run").disabled = false;
        }
      } catch (e) { /* 忽略网络抖动 */ }
    }, 2000);
  }

  function renderReport(job) {
    const report = (job.result && job.result.report_text) || "";
    const email = window.CONTACT_EMAIL || "xli@njucm.edu.cn";
    $("#result-body").innerHTML = `
      <div class="card" style="margin:0 0 14px;padding:14px 18px;background:#fefce8;border-color:#fde68a;">
        <p style="margin:0;font-size:13px;color:#92400e;">本平台仅提供文字报告。如需原始数据文件（复合物 PDB、轨迹等），请联系 <b>李雄博士（${email}）</b>。</p>
      </div>
      <div style="background:#f8fafc;border:1px solid var(--line);border-radius:10px;padding:18px 22px;font-size:14px;line-height:1.8;white-space:pre-wrap;font-family:'Segoe UI',system-ui,'Microsoft YaHei',sans-serif;">${esc(report)}</div>
    `;
    $("#result-card").classList.remove("result-hidden");
    $("#result-card").scrollIntoView({ behavior: "smooth" });
  }

  function esc(s) {
    return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  }

  document.addEventListener("DOMContentLoaded", init);
})();
