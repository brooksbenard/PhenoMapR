# Push scIMPEL to GitHub

Your project is ready to install from GitHub. To create the remote repo and push:

## 1. Create the repository on GitHub

1. Go to [github.com/new](https://github.com/new).
2. Set **Repository name** to `scIMPEL`.
3. Choose **Public**.
4. Do **not** add a README, .gitignore, or license (they already exist locally).
5. Click **Create repository**.

## 2. Add remote and push

In Terminal, from the `scIMPEL` directory:

```bash
# Replace YOUR_USERNAME with your GitHub username (e.g. bbenard)
git remote add origin https://github.com/YOUR_USERNAME/scIMPEL.git
git branch -M main
git push -u origin main
```

If you use SSH:

```bash
git remote add origin git@github.com:YOUR_USERNAME/scIMPEL.git
git push -u origin main
```

## 3. Install from GitHub (for you and others)

After pushing, anyone can install with:

```r
remotes::install_github("YOUR_USERNAME/scIMPEL")
```

Update the install command in `README.md` if your GitHub username is not `bbenard`.
