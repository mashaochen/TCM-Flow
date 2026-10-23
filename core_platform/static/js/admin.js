// DockLab 管理后台逻辑
(function () {
  "use strict";
  const $ = (s) => document.querySelector(s);
  window.logout = async function () { await fetch("/api/auth/logout", { method: "POST" }); location.href = "/"; };

  function esc(s) { return String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c])); }
  function notice(id, msg, type) { const el = $(id); el.className = "notice " + (type || "error"); el.textContent = msg; }

  async function loadStats() {
    const r = await fetch("/api/admin/stats"); const d = await r.json();
    if (!d.success) return;
    const s = d.stats;
    $("#stats").innerHTML = `
      <div class="stat"><div class="k">注册用户</div><div class="v">${s.n_users}</div></div>
      <div class="stat"><div class="k">对接任务</div><div class="v">${s.n_jobs}</div></div>
      <div class="stat"><div class="k">已完成</div><div class="v">${s.n_done}</div></div>
      <div class="stat"><div class="k">邀请码总数</div><div class="v">${s.n_codes}</div></div>
      <div class="stat"><div class="k">未使用邀请码</div><div class="v">${s.n_unused}</div></div>
    `;
  }

  async function loadInvites() {
    const r = await fetch("/api/admin/invites"); const d = await r.json();
    if (!d.success) return;
    const codes = d.codes;
    if (!codes.length) { $("#code-list").innerHTML = '<p style="color:#9ca3af;">暂无邀请码</p>'; return; }
    $("#code-list").innerHTML = codes.map((c) => `
      <span class="code-pill">
        <span class="mono">${esc(c.code)}</span>
        <span class="status ${esc(c.status)}">${statusText(c.status)}</span>
        ${c.used_by ? `<span style="font-size:11px;color:#6b7280;">${esc(c.used_by)}</span>` : ""}
        ${c.note ? `<span style="font-size:11px;color:#6b7280;">· ${esc(c.note)}</span>` : ""}
        ${c.expires_at ? `<span style="font-size:11px;color:#b45309;">· 至 ${esc(c.expires_at.slice(0, 10))}</span>` : `<span style="font-size:11px;color:#6b7280;">· 永久</span>`}
        <span style="font-size:11px;color:#9ca3af;">${esc(c.created_at)}</span>
        ${c.status === "unused" ? `<button style="border:none;background:none;cursor:pointer;color:#b91c1c;" onclick="revoke('${esc(c.code)}')">撤销</button>` : ""}
      </span>
    `).join("");
  }

  function statusText(s) {
    return s === "unused" ? "未使用" : s === "used" ? "已使用" : s === "revoked" ? "已撤销" : "已过期";
  }

  window.revoke = async function (code) {
    if (!confirm(`确定撤销邀请码 ${code}？`)) return;
    await fetch("/api/admin/invites/revoke", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ code }) });
    loadInvites(); loadStats();
  };

  async function loadUsers() {
    const r = await fetch("/api/admin/users"); const d = await r.json();
    if (!d.success) return;
    const users = d.users;
    if (!users.length) { $("#user-list").innerHTML = '<p style="color:#9ca3af;">暂无用户</p>'; return; }
    $("#user-list").innerHTML = `<table>
      <thead><tr><th>ID</th><th>用户名</th><th>角色</th><th>邀请码</th><th>注册时间</th></tr></thead>
      <tbody>${users.map((u) => `<tr><td>${u.id}</td><td>${esc(u.username)}</td>
        <td>${u.role === "admin" ? "管理员" : "用户"}</td>
        <td class="mono">${esc(u.invite_code || "-")}</td><td>${esc(u.created_at)}</td></tr>`).join("")}</tbody>
    </table>`;
  }

  async function loadJobs() {
    const r = await fetch("/api/admin/jobs"); const d = await r.json();
    if (!d.success) return;
    const jobs = d.jobs;
    if (!jobs.length) { $("#job-list").innerHTML = '<p style="color:#9ca3af;">暂无对接任务</p>'; return; }
    $("#job-list").innerHTML = `<table>
      <thead><tr><th>ID</th><th>用户</th><th>配体 SMILES</th><th>蛋白</th><th>状态</th><th>最佳亲和力</th><th>时间</th></tr></thead>
      <tbody>${jobs.map((j) => `<tr>
        <td class="mono">${j.id}</td><td>${esc(j.username)}</td>
        <td class="mono" style="max-width:180px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${esc(j.smiles)}">${esc(j.smiles)}</td>
        <td>${esc(j.pdb_id || "-")}</td>
        <td>${statusBadge(j.status)}</td>
        <td>${j.best_affinity != null ? j.best_affinity.toFixed(2) : "-"}</td>
        <td style="font-size:12px;">${esc(j.created_at)}</td></tr>`).join("")}</tbody>
    </table>`;
  }

  function statusBadge(s) {
    const map = { queued: "排队中", running: "运行中", done: "已完成", error: "失败" };
    return `<span class="status" style="background:#f1f5f9;color:#6b7280;padding:1px 8px;border-radius:10px;font-size:12px;">${map[s] || s}</span>`;
  }

  async function loadPackages() {
    const r = await fetch("/api/admin/packages"); const d = await r.json();
    if (!d.success) { $("#package-list").innerHTML = '<p style="color:#9ca3af;">加载失败</p>'; return; }
    const pkgs = d.packages || [];
    if (!pkgs.length) { $("#package-list").innerHTML = '<p style="color:#9ca3af;">暂无数据包（完成对接后会自动归档）</p>'; return; }
    let html = "";
    for (const day of pkgs) {
      html += `<h3 style="font-size:14px;margin:16px 0 6px;color:#0f766e;">📅 ${esc(day.date)}</h3>`;
      for (const u of day.users) {
        html += `<h4 style="font-size:13px;margin:8px 0 4px;color:#2563eb;">👤 ${esc(u.username)}（${u.jobs.length} 个任务）</h4>`;
        html += `<table><thead><tr><th>任务ID</th><th>配体</th><th>蛋白</th><th>亲和力</th><th>状态</th><th>数据文件（点击下载）</th></tr></thead><tbody>`;
        for (const j of u.jobs) {
          const dlLinks = (j.files || []).map(f =>
            `<a href="/api/admin/download/${j.job_id}/${f}" target="_blank" style="margin-right:6px;font-size:12px;">⬇ ${f}</a>`).join("");
          html += `<tr>
            <td class="mono">${j.job_id}</td>
            <td class="mono" style="max-width:150px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;" title="${esc(j.smiles)}">${esc(j.ligand_name || j.smiles)}</td>
            <td>${esc(j.pdb_id || "-")}</td>
            <td>${j.best_affinity != null ? j.best_affinity.toFixed(2) : "-"}</td>
            <td>${statusBadge(j.status)}</td>
            <td style="font-size:12px;">${dlLinks || '<span style="color:#9ca3af;">无文件</span>'}</td>
          </tr>`;
        }
        html += `</tbody></table>`;
      }
    }
    $("#package-list").innerHTML = html;
  }

  $("#btn-gen").addEventListener("click", async () => {
    notice("#gen-notice", "", "");
    const count = parseInt($("#invite-count").value) || 1;
    const note = $("#invite-note").value.trim();
    const validity = $("#invite-validity").value;
    const r = await fetch("/api/admin/invites", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ count, note, validity_days: validity === "" ? null : parseInt(validity) }),
    });
    const d = await r.json();
    if (!d.success) { notice("#gen-notice", d.error); return; }
    notice("#gen-notice", `成功生成 ${d.codes.length} 个邀请码`, "ok");
    $("#gen-result").innerHTML = d.codes.map((c) =>
      `<span class="code-pill"><span class="mono">${esc(c)}</span><span class="status unused">未使用</span></span>`
    ).join("");
    loadInvites(); loadStats();
  });

  async function init() {
    await Promise.all([loadStats(), loadInvites(), loadUsers(), loadJobs(), loadPackages()]);
  }
  document.addEventListener("DOMContentLoaded", init);
})();
