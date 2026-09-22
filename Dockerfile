# =============================================================================
# Dockerfile — eqtl-source-discordance-audit 可复现计算环境
# =============================================================================
# 用途: 从已处理的中游数据复现全部下游分析（方向一致性、Spearman 相关、
#       富集分析、Mahalanobis 匹配验证、跨人群复制评估、论文图表）
# =============================================================================
# 构建与运行:
#   docker build -t twas-eqtl-repro .
#   docker run --rm -v $(pwd)/output:/app/output twas-eqtl-repro
# =============================================================================

FROM continuumio/miniconda3:latest AS base

LABEL maintainer="Yijing Wu <1045381056@qq.com>"
LABEL description="eQTL-source discordance in TWAS — reproducibility environment (eqtl-source-discordance-audit)"
LABEL version="1.0"

# =============================================================================
# Stage 1: 系统依赖 + Conda 环境
# =============================================================================

# 安装系统级依赖（R 4.5.2 编译所需 + 实用工具）
RUN apt-get update && apt-get install -y --no-install-recommends \
    # R 基础环境
    r-base-core=4.5.2-1 \
    r-base-dev=4.5.2-1 \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    # 实用工具
    wget \
    ca-certificates \
    git \
    && rm -rf /var/lib/apt/lists/*

# 创建 conda 环境（Python 3.13 + 全部依赖）
COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && \
    conda clean -afy && \
    rm /tmp/environment.yml

# =============================================================================
# Stage 2: R 环境（matchit 包）
# =============================================================================

# 安装 matchit 及依赖的 R 包
RUN R -e "install.packages(c('MatchIt', 'cobalt', 'optmatch', 'Rcpp', 'RcppArmadillo'), \
    repos='https://cloud.r-project.org', \
    Ncpus=4, \
    clean=TRUE)" && \
    # 固定版本并保存 lockfile
    R -e "if(!require('renv')) install.packages('renv', repos='https://cloud.r-project.org'); \
    renv::init(bare=TRUE); \
    renv::install(c('MatchIt', 'cobalt', 'optmatch'), repos='https://cloud.r-project.org'); \
    renv::snapshot(type='all', lockfile='/app/renv.lock')"

# =============================================================================
# Stage 3: MetaXcan（S-PrediXcan 工具，用于完整管道复现）
# =============================================================================

RUN git clone --branch v0.8.1 https://github.com/hakyimlab/MetaXcan.git /opt/MetaXcan && \
    cd /opt/MetaXcan && \
    git checkout tags/v0.8.1 && \
    find /opt/MetaXcan -name "*.pyc" -delete

# 设置 MetaXcan 到 PATH
ENV PATH="/opt/MetaXcan/software:${PATH}"

# =============================================================================
# Stage 4: 项目文件
# =============================================================================

WORKDIR /app

# 复制项目文件（排除大文件）
COPY . /app/

# 创建输出目录（用于容器挂载）
RUN mkdir -p /app/output /app/data/raw /app/input

# 确保执行权限
RUN chmod +x /app/run_all.sh /app/run_spredixcan.sh 2>/dev/null || true

# =============================================================================
# 环境变量
# =============================================================================

# Conda 环境激活
ENV PATH="/opt/conda/envs/twas-eqtl/bin:${PATH}"
ENV CONDA_DEFAULT_ENV="twas-eqtl"

# 默认命令
ENTRYPOINT ["/bin/bash", "/app/run_all.sh"]
