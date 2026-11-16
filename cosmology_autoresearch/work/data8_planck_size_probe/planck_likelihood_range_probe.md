# Planck likelihood archive one-byte Range probe

Checked 2026-09-25 03:29 UTC. Issued exactly one HTTPS `GET` to the official ESA PLA product-action endpoint with `Range: bytes=0-0`. The server returned HTTP 200 (not 206), with no `Content-Range` or `Content-Length`, and `Transfer-Encoding: chunked`; it ignored the Range request. No redirect was returned, so the final URL is the requested URL. I stopped at the headers and read **zero response-body bytes**; no archive size can be inferred. Per the instruction, no retry or alternative request was made.

Response metadata:

- URL: `https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Data-baseline_R3.00.tar.gz`
- Status: `HTTP/1.1 200`
- `Content-Range`: absent
- `Content-Length`: absent
- `Transfer-Encoding`: `chunked`
- Redirects: none
- Final URL: same as requested
- Request bytes sent: 229 (application HTTP request)
- Response header bytes read: 291
- Response body bytes read: 0
- Archive size: unknown

This is not a download or package install. See [`planck_likelihood_range_probe.json`](planck_likelihood_range_probe.json) for machine-readable metadata.
