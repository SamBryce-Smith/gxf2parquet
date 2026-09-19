# Releasing gxf2parquet

This project publishes to [PyPI](https://pypi.org/project/gxf2parquet/)
automatically using [release-please](https://github.com/googleapis/release-please)
and PyPI [Trusted Publishing](https://docs.pypi.org/trusted-publishers/) (OIDC).

There are **two CI workflows**:

- `.github/workflows/ci.yml` — runs lint (`prek`/ruff) and the pytest matrix
  (Python 3.12 & 3.13) on every pull request and push to `main`. These are the
  checks `main` should require.
- `.github/workflows/release-please.yml` — on every push to `main`, maintains a
  "Release PR". Merging that PR bumps the version, updates `CHANGELOG.md`, tags
  the commit, creates a GitHub Release, and (in the same workflow) builds and
  uploads the package to PyPI.

## One-time setup

Do these once before the first release. **Step 1 must be done before the first
Release PR is merged**, or the publish step will fail.

### 1. Register the PyPI Trusted Publisher

This lets GitHub Actions upload to PyPI with no API token.

1. Sign in at <https://pypi.org>.
2. Go to **Your projects → Publishing** (for a brand-new project that does not
   exist on PyPI yet, use **Your account → Publishing → Add a pending publisher**:
   <https://pypi.org/manage/account/publishing/>).
3. Add a publisher with these **exact** values:
   - **PyPI Project Name:** `gxf2parquet`
   - **Owner:** `SamBryce-Smith`
   - **Repository name:** `gxf2parquet`
   - **Workflow filename:** `release-please.yml`
   - **Environment name:** `pypi`
4. Save.

> The owner, repo, workflow filename, and environment name must match the
> repository and `release-please.yml` exactly, or PyPI will reject the upload.

### 2. Create the `pypi` GitHub Environment

1. In the GitHub repo: **Settings → Environments → New environment**.
2. Name it exactly `pypi`.
3. Leave **required reviewers** empty — publishing is fully automatic.

### 3. Allow Actions to open the Release PR

1. **Settings → Actions → General**.
2. Under **Workflow permissions**, enable
   **"Allow GitHub Actions to create and approve pull requests"**.
3. Save.

### 4. Protect the `main` branch

1. **Settings → Rules → Rulesets → New branch ruleset** (or
   **Settings → Branches → Add branch protection rule**).
2. Target branch: `main`.
3. Enable:
   - **Require a pull request before merging.**
   - **Require status checks to pass** — add: `lint`, `test (3.12)`,
     `test (3.13)`.
   - **Require branches to be up to date before merging.**
4. Save.

> The status check names appear in the list only after the `ci.yml` workflow has
> run at least once (e.g. open this setup PR). If you don't see them yet, merge
> CI to `main` once, then add the checks.

## Day-to-day workflow

1. Open PRs against `main` using
   [Conventional Commit](https://www.conventionalcommits.org/) titles/messages:
   - `feat: ...` → minor version bump
   - `fix: ...` → patch version bump
   - `feat!: ...` or a `BREAKING CHANGE:` footer → major version bump
   - `chore: ...`, `docs: ...`, `ci: ...` → no release on their own
2. Merge your PR into `main` (CI must be green).
3. release-please opens or updates a **"chore(main): release x.y.z"** PR that
   stages the version bump and changelog.
4. When you're ready to ship, **merge the Release PR**. This tags the release,
   publishes the GitHub Release, and uploads to PyPI automatically.

## Verifying the first release

- After merging the Release PR, watch the **Release** workflow run; the `publish`
  job should appear and finish green.
- Confirm the new version at <https://pypi.org/project/gxf2parquet/>.
- Sanity check: `pip install gxf2parquet==<version>`.

## Notes / optional improvements

- **CI on the Release PR:** PRs opened by release-please use the built-in
  `GITHUB_TOKEN`, so their `pull_request` checks don't run automatically.
  Validation still happens via the `push` trigger on `main` after merge. If you
  want CI to run on the Release PR itself, pass a Personal Access Token or GitHub
  App token to the release-please action.
- **Enforce conventional commits on PR titles** (useful with squash-merge): add a
  PR-title linter such as
  [`amannn/action-semantic-pull-request`](https://github.com/amannn/action-semantic-pull-request).
