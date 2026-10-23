# -*- coding: utf-8 -*-
"""
Ollama API 兼容服务 — 用 llama-cpp-python 本地运行 Qwen2.5 模型，
提供与 Ollama 完全一致的 HTTP 接口，使 TCM-Flow 无需任何修改即可使用。

接口：
  GET  /api/tags        — 列出可用模型
  POST /api/generate    — 生成文本（支持 stream=False）
  POST /api/chat        — 聊天补全（可选）
"""
import os
import sys
import json
import time
from flask import Flask, request, jsonify

# ---------------------------------------------------------------------------
# 配置
# ---------------------------------------------------------------------------
MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models")
MODEL_NAME = os.environ.get("OLLAMA_MODEL_NAME", "qwen2.5:1.5b-instruct")
MODEL_FILE_PREFIX = os.environ.get("MODEL_FILE_PREFIX", "qwen2.5-1.5b-instruct-q4_k_m")
PORT = int(os.environ.get("OLLAMA_PORT", "11434"))
N_CTX = int(os.environ.get("N_CTX", "8192"))
N_GPU_LAYERS = int(os.environ.get("N_GPU_LAYERS", "0"))

# ---------------------------------------------------------------------------
# 加载模型
# ---------------------------------------------------------------------------
print("=" * 60)
print("  Ollama 兼容服务 · llama-cpp-python")
print("=" * 60)
print(f"模型目录: {MODEL_DIR}")
print(f"模型标识: {MODEL_NAME}")

model_path = None
for f in sorted(os.listdir(MODEL_DIR)):
    if f.startswith(MODEL_FILE_PREFIX) and f.endswith(".gguf"):
        model_path = os.path.join(MODEL_DIR, f)
        print(f"找到模型: {f}")
        break

if not model_path:
    print(f"错误: 未找到模型文件 (前缀: {MODEL_FILE_PREFIX})")
    sys.exit(1)

print("正在加载模型...")
t0 = time.time()
from llama_cpp import Llama

llm = Llama(
    model_path=model_path,
    n_ctx=N_CTX,
    n_gpu_layers=N_GPU_LAYERS,
    verbose=False,
)
print(f"模型加载完成，耗时 {time.time() - t0:.1f} 秒")

# 探测模型支持的 chat template（从模型元数据）
try:
    chat_template = llm.metadata.get("tokenizer.chat_template", "")
    print(f"Chat template: {chat_template[:80]}..." if len(chat_template) > 80 else f"Chat template: {chat_template}")
except Exception:
    print("Chat template: 未检测到，使用默认")

# ---------------------------------------------------------------------------
# Flask 服务
# ---------------------------------------------------------------------------
app = Flask(__name__)

SYSTEM_PROMPT = "你是一名专业的计算化学与分子模拟研究助理，擅长解读分子计算结果。回答要严谨、有据，不编造事实。对于推测性内容，明确标注为推测。"


def _generate_with_chat_template(user_prompt, system_prompt=None, temperature=0.3, max_tokens=2048):
    """
    使用 llama-cpp-python 内置的 create_chat_completion，自动应用模型的 chat template。
    这比手动组装更可靠，因为模板直接来自模型元数据。
    """
    messages = []
    if system_prompt:
        messages.append({"role": "system", "content": system_prompt})
    messages.append({"role": "user", "content": user_prompt})

    output = llm.create_chat_completion(
        messages=messages,
        max_tokens=max_tokens,
        temperature=temperature,
        # 让模型自动处理停止词
    )

    response_text = output["choices"][0]["message"]["content"].strip()
    usage = output.get("usage", {})
    return response_text, usage


@app.route("/api/tags", methods=["GET"])
def api_tags():
    return jsonify({
        "models": [{
            "name": MODEL_NAME,
            "model": MODEL_NAME,
            "size": os.path.getsize(model_path),
            "modified_at": "2024-01-01T00:00:00Z",
            "details": {
                "format": "gguf",
                "family": "qwen2",
                "parameter_size": "1.5B",
                "quantization_level": "Q4_K_M",
            }
        }]
    })


@app.route("/api/generate", methods=["POST"])
def api_generate():
    """
    Ollama 兼容的 /api/generate 接口。
    TCM-Flow llm_report.py 通过此接口调用。
    """
    data = request.get_json()
    prompt = data.get("prompt", "")
    stream = data.get("stream", False)
    options = data.get("options", {})
    temperature = options.get("temperature", 0.3)
    num_predict = options.get("num_predict", 2048)
    
    t0 = time.time()
    response_text, usage = _generate_with_chat_template(
        prompt,
        system_prompt=SYSTEM_PROMPT,
        temperature=temperature,
        max_tokens=num_predict,
    )
    elapsed = time.time() - t0
    
    eval_count = usage.get("completion_tokens", 0)
    
    result = {
        "model": MODEL_NAME,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "response": response_text,
        "done": True,
        "total_duration": int(elapsed * 1e9),
        "load_duration": 0,
        "prompt_eval_count": usage.get("prompt_tokens", 0),
        "prompt_eval_duration": 0,
        "eval_count": eval_count,
        "eval_duration": int(elapsed * 1e9),
    }
    
    return jsonify(result)


@app.route("/api/chat", methods=["POST"])
def api_chat():
    """可选：Ollama 兼容的 /api/chat 接口。"""
    data = request.get_json()
    messages = data.get("messages", [])
    stream = data.get("stream", False)
    options = data.get("options", {})
    temperature = options.get("temperature", 0.3)
    num_predict = options.get("num_predict", 2048)

    t0 = time.time()
    output = llm.create_chat_completion(
        messages=messages,
        max_tokens=num_predict,
        temperature=temperature,
    )
    elapsed = time.time() - t0

    response_text = output["choices"][0]["message"]["content"].strip()
    usage = output.get("usage", {})

    result = {
        "model": MODEL_NAME,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "message": {
            "role": "assistant",
            "content": response_text,
        },
        "done": True,
        "total_duration": int(elapsed * 1e9),
        "load_duration": 0,
        "prompt_eval_count": usage.get("prompt_tokens", 0),
        "prompt_eval_duration": 0,
        "eval_count": usage.get("completion_tokens", 0),
        "eval_duration": int(elapsed * 1e9),
    }

    return jsonify(result)


@app.route("/", methods=["GET"])
def index():
    return jsonify({
        "service": "Ollama-compatible LLM server",
        "model": MODEL_NAME,
        "backend": "llama-cpp-python",
        "status": "running",
    })


if __name__ == "__main__":
    print(f"启动服务: http://localhost:{PORT}")
    print("=" * 60)
    app.run(host="0.0.0.0", port=PORT, debug=False)
