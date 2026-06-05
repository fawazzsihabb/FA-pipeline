# Publishing this repo to GitHub

I can't push to your GitHub account for you (no credentials), so here are the two reliable ways. The folder is already structured as a complete repo (README, src/, notebooks/, data/, figures/, LICENSE, .gitignore).

## Option A — GitHub CLI (fastest, keeps folder structure)
Install gh once: https://cli.github.com  (macOS: `brew install gh`, Ubuntu: `sudo apt install gh`)
From inside the unzipped folder:
```bash
gh auth login
git init
git add .
git commit -m "FA19 COX-2 computational drug-discovery portfolio"
gh repo create ferulic-cox2-screening --public --source=. --remote=origin --push
```
Done — it prints your repo URL.

## Option B — plain git with an existing empty repo
Create an empty repo on github.com first (no README), then:
```bash
git init
git add .
git commit -m "FA19 COX-2 computational drug-discovery portfolio"
git branch -M main
git remote add origin https://github.com/<your-username>/ferulic-cox2-screening.git
git push -u origin main
```

## Option C — GitHub Desktop (no terminal)
Open GitHub Desktop → File → Add Local Repository → pick the unzipped folder →
"Publish repository". (The website's drag-and-drop upload flattens folders, so prefer Desktop or the CLI.)

Notebooks with embedded figures are fine (well under GitHub's 100 MB/file limit).
