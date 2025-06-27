# DGscRNA Package Publishing Guide

This guide explains how to publish the DGscRNA package to PyPI using the automated GitHub workflows.

## 🚀 Quick Start

1. **Set up PyPI accounts and tokens** (one-time setup)
2. **Configure GitHub repository secrets** (one-time setup)
3. **Test on TestPyPI** (recommended before first release)
4. **Create a release** to automatically publish to PyPI

## 📋 Prerequisites

### PyPI Account Setup

1. **Create accounts**:
   - Register at [PyPI](https://pypi.org/account/register/)
   - Register at [TestPyPI](https://test.pypi.org/account/register/) (for testing)

2. **Generate API tokens**:
   - **PyPI**: Go to [PyPI Account Settings](https://pypi.org/manage/account/) → API tokens → Add API token
   - **TestPyPI**: Go to [TestPyPI Account Settings](https://test.pypi.org/manage/account/) → API tokens → Add API token
   - Save these tokens securely (they start with `pypi-`)

### GitHub Repository Setup

1. Go to your GitHub repository
2. Navigate to **Settings** → **Secrets and variables** → **Actions**
3. Add these repository secrets:

| Secret Name | Value | Description |
|-------------|-------|-------------|
| `PYPI_API_TOKEN` | Your PyPI token | For publishing to PyPI |
| `TEST_PYPI_API_TOKEN` | Your TestPyPI token | For testing on TestPyPI |

## 🔄 Workflows Overview

### 1. CI Workflow (`.github/workflows/ci.yml`)
- **Triggers**: Push/PR to main/develop branches
- **What it does**:
  - Tests across Python 3.8-3.11 on Ubuntu/Windows/macOS
  - Runs linting with flake8
  - Tests package build and installation
  - Runs pytest (basic import tests)

### 2. Publish Workflow (`.github/workflows/python-publish.yml`)
- **Triggers**: 
  - Automatic: When a new release is published
  - Manual: Via workflow dispatch
- **What it does**:
  - Builds the package (wheel and source distribution)
  - Tests installation across multiple platforms
  - Publishes to TestPyPI or PyPI

## 📦 Publishing Process

### Step 1: Test Build Locally (Optional but Recommended)

```bash
# Install build tools
pip install build twine

# Build the package
python -m build

# Check the built package
twine check dist/*

# Test installation locally
pip install dist/*.whl
python -c "import dgscrna; print('Success!')"
```

### Step 2: Test on TestPyPI

1. Go to **Actions** tab in your GitHub repository
2. Click on **"Publish to PyPI"** workflow
3. Click **"Run workflow"**
4. Select:
   - **Use workflow from**: `main`
   - **Environment to deploy to**: `testpypi`
5. Click **"Run workflow"**

This will build and upload to TestPyPI. You can then test installation:

```bash
pip install -i https://test.pypi.org/simple/ dgscrna
```

### Step 3: Publish to PyPI

#### Method 1: Automatic (Recommended)

1. **Create a new release**:
   - Go to your repository's main page
   - Click **"Releases"** → **"Create a new release"**
   - Click **"Choose a tag"** → Type new version (e.g., `v0.1.1`)
   - **Release title**: `Release v0.1.1`
   - **Description**: Add release notes describing changes
   - Click **"Publish release"**

2. **Automatic workflow**:
   - The release automatically triggers the publish workflow
   - Package builds and publishes to PyPI
   - Check the Actions tab for progress

#### Method 2: Manual

1. Go to **Actions** → **"Publish to PyPI"**
2. Click **"Run workflow"**
3. Select `pypi` environment
4. Click **"Run workflow"**

## 🔢 Version Management

### Update Version Number

Before publishing, update the version in `pyproject.toml`:

```toml
[project]
name = "dgscrna"
version = "0.1.1"  # Update this line
```

### Semantic Versioning

Follow [semantic versioning](https://semver.org/):
- **Patch** (0.1.0 → 0.1.1): Bug fixes, no new features
- **Minor** (0.1.0 → 0.2.0): New features, backward compatible
- **Major** (0.1.0 → 1.0.0): Breaking changes

## ✅ Pre-Publication Checklist

- [ ] Update version number in `pyproject.toml`
- [ ] Update `README.md` if needed
- [ ] Test package locally
- [ ] Test on TestPyPI first
- [ ] All CI tests pass
- [ ] Write meaningful release notes

## 🔍 Verification

After publishing, verify the publication:

1. **Check PyPI page**: https://pypi.org/project/dgscrna/
2. **Test installation**:
   ```bash
   pip install dgscrna
   python -c "import dgscrna; print('Published successfully!')"
   ```

## 🐛 Troubleshooting

### Common Issues

1. **"Package already exists"**:
   - The version was already published
   - Bump the version number and try again

2. **"Invalid token"**:
   - Check that tokens are correctly set in repository secrets
   - Regenerate tokens if needed

3. **"Build failed"**:
   - Check the CI workflow logs for specific errors
   - Ensure all dependencies are correctly specified
   - Fix any import or syntax errors

4. **"Tests failed"**:
   - Review test output in the Actions tab
   - Fix failing tests before publishing

### Debug Steps

1. **Check workflow logs**: Go to Actions tab → Click on failed workflow → Review logs
2. **Test locally**: Run the same commands locally to reproduce issues
3. **Check dependencies**: Ensure all required packages are in `pyproject.toml`
4. **Verify package structure**: Ensure all modules can be imported

## 📁 Package Structure

Your package should have this structure:
```
DGscRNA/
├── dgscrna/                 # Main package
│   ├── __init__.py
│   ├── core/               # Core modules
│   └── models/             # Model modules
├── tests/                  # Test files
├── pyproject.toml          # Modern packaging config
├── setup.py               # Backward compatibility
├── README.md              # Package description
├── LICENSE                # GPL-3.0 license
└── .github/workflows/     # GitHub Actions
```

## 🔒 Security Best Practices

- Never commit API tokens to the repository
- Use repository secrets for sensitive information
- Regularly rotate API tokens
- Monitor package downloads and usage
- Review dependencies for security vulnerabilities

## 📈 Post-Publication

After successful publication:

1. **Monitor**: Check download statistics on PyPI
2. **Maintain**: Keep dependencies updated
3. **Support**: Respond to user issues on GitHub
4. **Iterate**: Continue development and release new versions

## 🆘 Getting Help

If you encounter issues:

1. Check the [GitHub Actions documentation](https://docs.github.com/en/actions)
2. Review [PyPI publishing guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
3. Check the workflow logs for detailed error messages
4. Ensure package structure matches Python packaging standards

---

**Example Complete Workflow:**

```bash
# 1. Local development
git checkout main
git pull origin main

# 2. Update version
# Edit pyproject.toml: version = "0.1.1"

# 3. Test locally
python -m build
pip install dist/*.whl

# 4. Test on TestPyPI (via GitHub Actions)
# Actions → Publish to PyPI → Run workflow (testpypi)

# 5. Create release (automatic publishing)
# Releases → Create new release → v0.1.1

# 6. Verify
pip install dgscrna
python -c "import dgscrna; print('Success!')"
```

This automated workflow ensures your DGscRNA package is thoroughly tested and properly published to PyPI! 🎉 