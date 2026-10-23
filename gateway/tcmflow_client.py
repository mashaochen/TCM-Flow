# -*- coding: utf-8 -*-
"""
TCM-Flow API 客户端
通过HTTP调用TCM-Flow平台的全流程计算接口
TCM-Flow本身代码完全不修改
"""
import requests
import json
import config


class TCMFlowClient:
    """TCM-Flow 平台 API 客户端"""
    
    def __init__(self):
        self.base_url = config.TCMFLOW_BASE_URL.rstrip("/")
        self.session = requests.Session()
        self._logged_in = False
    
    def _ensure_login(self):
        """确保已登录"""
        if not self._logged_in:
            self._login()
    
    def _login(self):
        """登录TCM-Flow"""
        try:
            resp = self.session.post(
                f"{self.base_url}/api/auth/login",
                json={
                    "username": config.TCMFLOW_ADMIN_USER,
                    "password": config.TCMFLOW_ADMIN_PASS
                },
                timeout=10
            )
            data = resp.json()
            if data.get("success"):
                self._logged_in = True
                return True
            return False
        except Exception as e:
            print(f"[TCMFlowClient] 登录失败: {e}")
            return False
    
    def run_full_pipeline(self, smiles, ligand_name, pdb_id=None, 
                          pdb_text=None, dock_params=None, manual_box=None):
        """
        提交全流程计算任务
        
        参数:
            smiles: 配体SMILES
            ligand_name: 配体名称
            pdb_id: PDB ID（三选一）
            pdb_text: PDB文本内容（三选一，与上传文件等价）
            dock_params: 对接参数字典
            manual_box: 手动对接盒
        
        返回: (success, job_id/error_msg)
        """
        self._ensure_login()
        
        # 判断受体来源
        if pdb_id:
            receptor_source = "pdb_id"
        elif pdb_text:
            receptor_source = "pdb_text"
        else:
            return False, "请提供PDB ID或PDB文本"
        
        payload = {
            "smiles": smiles,
            "ligand_name": ligand_name,
            "receptor_source": receptor_source,
            "pdb_id": pdb_id or "",
            "pdb_text": pdb_text or "",
            "dock_params": dock_params or {},
        }
        if manual_box:
            payload["manual_box"] = manual_box
        
        try:
            resp = self.session.post(
                f"{self.base_url}/api/full/run",
                json=payload,
                timeout=30
            )
            data = resp.json()
            if data.get("success"):
                return True, data["job_id"]
            else:
                return False, data.get("error", "未知错误")
        except Exception as e:
            return False, f"调用TCM-Flow失败: {e}"
    
    def get_job_status(self, job_id):
        """查询任务状态"""
        self._ensure_login()
        
        try:
            resp = self.session.get(
                f"{self.base_url}/api/full/{job_id}",
                timeout=10
            )
            data = resp.json()
            if data.get("success"):
                return True, data["job"]
            else:
                return False, data.get("error", "未知错误")
        except Exception as e:
            return False, f"查询失败: {e}"
    
    def health_check(self):
        """检查TCM-Flow是否在线"""
        try:
            resp = requests.get(f"{self.base_url}/api/health", timeout=5)
            data = resp.json()
            return data.get("status") == "ok", data
        except Exception as e:
            return False, {"error": str(e)}


# 全局单例
_client = None

def get_client():
    global _client
    if _client is None:
        _client = TCMFlowClient()
    return _client
