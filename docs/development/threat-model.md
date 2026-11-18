# Threat model

bijli is a computation library: most of its API takes numbers the caller chose and returns numbers. Three surfaces
read data that may come from elsewhere, and one talks to the network. This file names them, what each defends
against, and what remains (M8 part 2).

## Surfaces

| Surface | Input | Reached by |
|---|---|---|
| The Touchstone reader | a Touchstone v1.1 file's bytes and length | `bj_parse_touchstone` |
| The JSON readers | a JSON document's bytes and length | every `BjX_from_json` (`src/json.cyr`, `src/ai.cyr`) |
| The soorat constructors | a component name's bytes and length | `BjFieldSlice2D_from_fdtd_2d` and the 3D pair |
| The ai client | daimon's HTTP answer; the configured endpoint and API key | `BjDaimonClient::register_agent` (`[lib.ai]` only) |

## What each defends against

- **Malformed text.** Every byte input is an address and a length, never NUL-terminated, and is checked to be UTF-8
  (no overlong forms, surrogates or values past U+10FFFF) before any of it is read as text (ADR 0002 rows 38 and 44).
  A refusal returns an error code with the reason and byte offset in the last-error record, and leaves `out`
  untouched.
- **Unbounded documents.** The Touchstone reader takes 1 to 256 ports (Rust's `MAX_TOUCHSTONE_PORTS`) and sizes its
  one block from the records it read. A JSON document is parsed once into a tree in a decode region (freelist chunks
  on the reader's stack), given back before the reader returns; nesting past 128 levels, a key repeated in any
  object, a lone surrogate escape or a raw control character anywhere refuses the whole document (row 12), and a
  count past 2⁶³ − 1 is refused (row 36). A block that can't be had is `BJ_ERR_ALLOCATION`, never an abort (row 9).
- **Numbers.** Every number is read correctly rounded (row 11; an f32 rounded once, row 45); one whose value
  overflows is refused. NaN and ±∞ can't be written as JSON numbers, so a writer emits `null` and the reader refuses
  it.
- **The ai client's requests.** An endpoint or API key holding a control byte (U+0000..U+001F, U+007F) is refused
  before anything is sent, so neither can split a request or inject a header. A request has a 30 s deadline and
  follows no redirect. The response is capped (a 1 MiB body, 1.06 MiB in all); a status outside 2xx, an oversized,
  malformed or truncated answer, and one without an integer `id` are all `BJ_ERR_TRANSPORT`. The exchange's memory
  comes from a decode region freed before the call returns. TLS (an `https` endpoint) is sandhi's, with its policy.
- **Secrets.** The API key is never logged (the debug lines print the endpoint only), and `debug` prints it as
  `Some("<redacted>")`. The JSON writer writes it in clear, as Rust's `Serialize` did: a written config is a secret.
- **Memory.** No time step allocates; a grid, matrix or slice owns its blocks and its `_free` gives them back once (a
  second call does nothing). The Touchstone reader is fuzzed (`tests/bijli.fcyr`: byte-level mutants of every format,
  1 to 4 ports, checked to return, to leave `out` untouched on a refusal, and not to grow the heap on a replay).

## What remains

- **The JSON readers have no fuzz harness.** They rest on bayan's parser and on the golden corpus's documents
  (malformed, repeated keys, deep nesting, every escape, numbers at every edge), not on mutation. A harness like
  `tests/bijli.fcyr` for `_bj_jr_doc` and a few readers would close it.
- **Size is the caller's to bound.** The readers take whatever length they are given; a multi-gigabyte document is
  parsed (or refused as Allocation) as asked. A service reading untrusted documents should cap their size first.
- **The ai client trusts its endpoint.** It authenticates daimon only as far as TLS does for an `https` endpoint; an
  `http` endpoint is plaintext, the API key included.
- **Numerical results are not a security boundary.** A wrong value is a bug to report as one, not a vulnerability,
  unless it reads or writes memory it shouldn't.
