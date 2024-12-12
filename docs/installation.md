# Pangenome Heritability Tool

A tool for processing pangenome structural variants and generating PLINK format files.

## Installation Guide

### Requirements
- Python 3.8+
- MUSCLE 5

### Step 1: Install MUSCLE
```bash
# Linux
wget https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3
chmod +x muscle-linux-x86.v5.
# Add your local bin to PATH
# Add these lines to your ~/.bashrc (for bash) or ~/.zshrc (for zsh)
echo 'export PATH="$HOME/local/bin:$PATH"' >> ~/.bashrc
# OR for zsh
# echo 'export PATH="$HOME/local/bin:$PATH"' >> ~/.zshrc

# Reload your shell configuration
source ~/.bashrc  # OR source ~/.zshrc for zsh

# macOS
brew install muscle
```

### Step 2: Install Python Package
```bash
# Clone repository
git clone https://github.com/yuanpeixiong/pangenome-heritability.git
cd pangenome-heritability

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/macOS
.\venv\Scripts\activate   # Windows

# Install dependencies
pip install -e .
```




