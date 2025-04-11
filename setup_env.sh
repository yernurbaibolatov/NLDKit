#!/bin/bash

echo "🔧 Setting up Python virtual environment for NLDKit..."

# Step 1: Create the virtual environment
python3 -m venv .venv
echo "✅ Virtual environment created at .venv/"

# Step 2: Activate the environment
source .venv/bin/activate

# Step 3: Upgrade pip
echo "⬆️  Upgrading pip..."
pip install --upgrade pip

# Step 4: Install runtime dependencies
echo "📦 Installing core dependencies..."
pip install numpy scipy matplotlib

# Step 5: Install development tools
echo "🛠️  Installing dev tools: pytest, black, mypy, jupyter"
pip install pytest black mypy jupyter ipykernel

# Step 6: Freeze environments
pip freeze > requirements.txt
pip freeze > requirements-dev.txt
echo "📄 requirements.txt and requirements-dev.txt created."

# Step 7: Configure VS Code Python interpreter
mkdir -p .vscode
cat > .vscode/settings.json << EOF
{
  "python.pythonPath": "\${workspaceFolder}/.venv/bin/python"
}
EOF
echo "🧠 VS Code interpreter set."

# Step 8: Register Jupyter kernel
python -m ipykernel install --user --name=nldkit --display-name "Python (NLDKit)"
echo "📚 Jupyter kernel 'Python (NLDKit)' installed."

echo "🎉 All done! Activate with: source .venv/bin/activate"
