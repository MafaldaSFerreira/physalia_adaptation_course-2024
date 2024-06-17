# Useful commands on Git/GitHub

## On Git and GitHub
For a brief introduction to Git/GitHub developed by Eric Normandeau, check this [webpage](https://github.com/enormandeau/github_tutorial).

For a more extended introduction, see this [webpage](https://docs.github.com/en/get-started/using-git/about-git#how-github-works).

If you use [Visual Studio Code](https://code.visualstudio.com) as your code editor, note that it has specific features that enable [using Git within VS Code](https://code.visualstudio.com/docs/sourcecontrol/intro-to-git).


## before starting
```bash
# to clone the repository the first time to your local computer
git clone
# to update the repository when you open it
git pull
```

## during analysis
If you add big files/directories to your local copy, update the `.gitignore` file to avoid adding them to the GitHub repository when syncing.

If you want to upload an empty folder, create the folder and add inside a `README.MD` file with one line explaining what that directory will contain.

## finishing editions
```bash
# to add the files that you have made or update
git add .
# to document what are the changes you did
git commit -m "blabla"
# to update the files you add to the web repository
git push
```
