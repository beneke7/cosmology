# Planck PR3 archive staging report

Status: the exact official Planck PR3 baseline likelihood archive was downloaded and validated. Its complete gzip stream and tar structure are readable. The archive was not extracted, and no packages or likelihoods were run.

## Transfer

- Requested and final URL: `https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_Likelihood_Data-baseline_R3.00.tar.gz`
- HTTP: 200; redirects: 0; content type: `application/x-gzip`; content disposition names `COM_Likelihood_Data-baseline_R3.00.tar.gz`.
- Transfer encoding: chunked; no `Content-Length` was sent. No cookie jar, account, or authentication credentials were configured.
- Exact archive size: 60,323,470 bytes. SHA-256: `0b73171e3acc671c28184466a45485a2d1c1d93676b832abdfe688c7b04024e6`.
- Command: `timeout 7200s python3 -u cosmology_autoresearch/work/planck_stage/download_capped.py` from `/home/v/proj/bene/cosmology`. The Python standard-library downloader streams in blocks and stops at the authorized 21,365,665,738-byte cap. Transfer time was 6.028 seconds.
- Response headers are in [transfer.headers.txt](transfer.headers.txt); the byte-count and hash log is in [transfer.log](transfer.log).

## Archive inventory

The gzip stream passed a full read through EOF. The uncompressed tar stream is 92,425,728 bytes. The tar contains 362 members: 294 regular files totaling 92,163,979 declared file bytes and 68 directories. No links or other member types were present. No files were extracted.

The three products identified in the Planck gate are present as `.clik` directory trees:

| Likelihood | Archive member root | Members | Regular-file bytes |
| --- | --- | ---: | ---: |
| Plik HM TTTEEE | `baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik` | 138 | 57,891,282 |
| Commander low-TT | `baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik` | 11 | 6,509,031 |
| SimAll lowE EE | `baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik` | 10 | 694,285 |

The complete tar member list is [tar_members.jsonl](tar_members.jsonl). The listing and per-likelihood member details, byte totals, transfer metadata, and cap arithmetic are also recorded in [inspection.json](inspection.json).

## Cache-cap arithmetic

The selected-cache limit is 21,474,836,480 bytes (20 GiB), including 109,170,742 existing manifest bytes. That leaves the authorized 21,365,665,738-byte maximum archive transfer. After staging this 60,323,470-byte archive, existing cache plus archive is 169,494,212 bytes, leaving 21,305,342,268 bytes.

Counting both the retained archive and the full regular-file payload as if extracted gives 261,658,191 bytes total, below the 20 GiB cap. Extraction was not part of this task and was not performed.

## Rights and limits

The archive is distributed under ESA's CC BY-NC 3.0 IGO terms. There were no transfer or format failures. The endpoint did not provide a content length; the streaming byte cap controlled that case. This report covers only archive staging and inspection, not package installation or likelihood execution.
