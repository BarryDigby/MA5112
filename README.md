# MA5112
2023 MA5112

## Github

You should already have the MA5112 git repository on your laptop from last friday. Change to the directory and pull the latest changes I have published:

```bash
cd MA5112/
git pull origin main
```

### Authentication tokens

You need to set up a token for CLI access to publish your branch to GitHub. This is why you were all getting errors during last fridays practical.

Follow the instructions at this page: [https://dev.to/ibmdeveloper/can-t-push-to-your-github-repo-i-can-help-with-that-1fda](https://dev.to/ibmdeveloper/can-t-push-to-your-github-repo-i-can-help-with-that-1fda).

**You do not need to use the Keychain access application in the tutorial above. Simply save the token somewhere as you will be asked to produce it when pushing to GitHub for the first time (along with your username and password). You will not be prompted for this information again. Keys can be regenerated.**

### Publish your branch

**This will not work until you have accepted my collaboration reqest for the MA5112 repository on GitHub.**

Run the below to stage your commits (your work last friday) and publish them to the MA5112 repo under your branch name.

```console
cd MA5112/
git branch
```

Make sure you are on your branch and not `main` (run `git checkout -b <your branch name>` if needed)

```console
git add . 
git commit -m "my first commit"
git push -u origin <your branch name>
```

### Troubleshooting

When running `git pull origin main` in the future you are likely to encounter this:

```console
From github.com:BarryDigby/MA5112
 * branch            main       -> FETCH_HEAD
hint: You have divergent branches and need to specify how to reconcile them.
hint: You can do so by running one of the following commands sometime before
hint: your next pull:
hint: 
hint:   git config pull.rebase false  # merge
hint:   git config pull.rebase true   # rebase
hint:   git config pull.ff only       # fast-forward only
hint: 
hint: You can replace "git config" with "git config --global" to set a default
hint: preference for all repositories. You can also pass --rebase, --no-rebase,
hint: or --ff-only on the command line to override the configured default per
hint: invocation.
fatal: Need to specify how to reconcile divergent branches.
```

We want to use `git config pull.rebase false` to accept the latest changes from `main` without overwriting your own work.