# Hosting the hub

A UCSC track hub must be reachable over **public HTTP/HTTPS with no authentication**.
UCSC fetches `hub.txt` and every file it references by relative path, so upload the whole
output directory intact.

## Hosting options

| Option | Notes |
|--------|-------|
| **GitHub (raw URL)** | Free, easy. Repo must be **public**. Load the RAW URL, not the repo page (see below). |
| **Institutional / lab server** | `scp -r my_hub/ user@server.edu:/public_html/hubs/` → `http://server.edu/~user/hubs/my_hub/hub.txt`. Must allow direct file access, no auth. |
| **AWS S3** | `aws s3 cp my_hub/ s3://my-bucket/my_hub/ --recursive --acl public-read` → `https://my-bucket.s3.amazonaws.com/my_hub/hub.txt`. |
| **Google Cloud Storage** | `gsutil -m cp -r my_hub/ gs://my-bucket/` then `gsutil -m acl ch -u AllUsers:R gs://my-bucket/my_hub/**` → `https://storage.googleapis.com/my-bucket/my_hub/hub.txt`. |

### GitHub: use the RAW URL

```
Correct: https://raw.githubusercontent.com/USER/REPO/BRANCH/hub.txt
Wrong:   https://github.com/USER/REPO/blob/BRANCH/hub.txt
```

Click `hub.txt` in the repo → "Raw" button → copy the URL, or construct it as above.

## Not recommended

- **Dropbox / Google Drive** — require auth, no direct file access.
- **Private repositories** — UCSC cannot read auth-protected files.

## Validate before uploading

```bash
hubCheck -noTracks hub.txt
```

`-noTracks` checks hub/genome/trackDb structure without downloading every bigBed — fast and
enough to catch path/config errors. You can also validate a remote hub after upload:
`hubCheck -noTracks https://.../hub.txt`. Validate inputs before you even build the hub with
`python -m sqanti_browser ... --validate-only`.

## Load in UCSC

1. Go to https://genome.ucsc.edu/
2. **My Data → Track Hubs**
3. Open the **Connected Hubs** tab and paste the hub URL (must point to `hub.txt`)
4. Click **Add Hub**
5. Select the genome and navigate to view tracks

First load may show tracks in "hide" mode — set them to "full" or "pack".

## Updating a hosted hub

After re-uploading regenerated files, UCSC caches aggressively. Disconnect the hub, add
`udcTimeout=5` to the browser URL (`.../hgTracks?udcTimeout=5`) to force a cache refresh,
then reconnect.

## Common issues

| Problem | Fix |
|---------|-----|
| Hub won't load | Confirm URL is public (try incognito); `hub.txt` at the URL you pasted; run `hubCheck -noTracks hub.txt`. |
| 404 on `.bb` files | Ensure every referenced file was uploaded; trackDb uses relative paths. |
| Changes not appearing | Disconnect, add `udcTimeout=5`, reconnect. |
| GitHub wrong format | Must be `raw.githubusercontent.com`. |
