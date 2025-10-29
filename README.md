# BSProject2526
💻 Repository Workflow Guide

This guide outlines the standard Git workflow for collaborating on this repository, focusing on using a Personal Access Token (PAT) for authentication (required by GitHub over HTTPS) and the essential three-step cycle: Stage, Commit, Push.

1. Authentication Setup (Using a Personal Access Token)

GitHub requires a Personal Access Token (PAT) instead of your regular password for all command-line operations (like git push) over HTTPS.

🔑 A. Generate the PAT

Go to GitHub: Navigate to Settings > Developer settings > Personal access tokens > Tokens (classic).

Generate a New Token: Click Generate new token (classic).

Configure:

Note: Give it a meaningful name (e.g., My-Work-PC).

Expiration: Set an appropriate expiration date (e.g., 90 days or 1 year).

Scope: You must check the repo box. This grants the token read/write access to your repositories.

Save the Token: Once generated, IMMEDIATELY COPY the resulting string. This is the only time you will ever see it.

💾 B. Use and Store the PAT

The first time you run a command that requires authentication (git push or git clone for a private repo), Git will prompt you for credentials:

Username: Enter your GitHub username.

Password: PASTE THE PAT (the long string) here, NOT your actual GitHub account password.

Your system's Git Credential Manager will securely store this PAT so you will not be prompted for it again until it expires.

2. Initial Setup: Cloning the Repository

Cloning downloads a copy of the repository to your local computer and sets up the necessary links to GitHub.

Navigate: Open your terminal or Git Bash and move to the directory where you want to store the project:

cd ~/Projects/


Clone: Copy the HTTPS or SSH link from the GitHub repository page and clone it:

git clone [repository_url]


Enter Directory:

cd [repository-name]


3. Daily Workflow: Making and Pushing Changes

The core Git workflow is a three-step cycle that turns local file modifications into recorded history on GitHub.

Step 1: Check Status and Modify Files

Before you start, ensure your local copy is up-to-date and check what you're working on.

Command

Purpose

git pull

Downloads the latest changes from GitHub to ensure you have the most recent version.

git status

Shows which files are modified (edited), untracked (new), or staged (ready to commit).

Action

Create, Edit, or Delete the files in your repository folder.

Step  Step 2: Stage and Commit (Saving the Change)

You must tell Git which changes you want to include in the next version snapshot.

📁 A. Adding / Editing a File

To include new files or modified files in your commit, you must Stage them first:

# To stage a single file:
git add my_new_file.txt

# To stage ALL modified, new, and deleted files:
git add .


✍️ B. Committing the Snapshot

Once files are staged, you create a permanent snapshot (a commit) in your local history:

git commit -m "Brief but descriptive summary of the changes"


🗑️ C. Deleting a File

If you want to remove a file from the repository, you must use git rm to stage the deletion:

# Deletes the file locally AND stages the deletion for the next commit
git rm file_to_delete.txt 

# Then commit the deletion:
git commit -m "Removed old file_to_delete.txt"


Step 3: Push to Remote (Updating GitHub)

This command uploads all the new commits from your local repository to GitHub.

git push


If this is the first time you've pushed a new branch, you may need to use:

git push -u origin [branch-name]


⚠️ Troubleshooting Tip

If your file is present locally but not on GitHub, you likely missed the git add step. Always run git status to verify your file is in the "Changes to be committed" (green) section before committing.
