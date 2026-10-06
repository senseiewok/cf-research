# Review: YouTube channel About text against the lab's own rules

Status: proposal. Channel `@SenseiEwok`, About page read 2026-10-04. No videos yet, four external links (the lab website, senseiewok.ai; CFF WA donation page, CFF WA chapter page, cff.org). This is the only public text the channel has, so it is worth getting right before the first upload.

The review standard is not my taste. It is `hq/README.md`, `cf-research-context`, and `lab-voice`, applied to the channel as they would be to a README.

## Where the current text conflicts with the repos

| Current About text | Rule it conflicts with | Why it matters |
| --- | --- | --- |
| "accelerate the cure for Cystic Fibrosis" | `hq/README.md`: "We do not make medicine, provide clinical care, or claim a cure." `lab-voice`: "Never promise a cure, access to care, a timeline, or a clinical benefit from our experiments." | This is the single sentence a journalist, a clinician, or the Foundation would quote back. The repos say the opposite in plain words |
| "our medical AI projects" | `cf-research-context`: the lab builds research tooling and "may not ... present a tool as a replacement for clinical care" | "Medical AI" is read as clinical software. The repos describe schemas, literature tools, and claim-checking |
| "designed to ease the burden on the CF community" | `lab-voice`: "Hope without a forecast"; separate aspirations from tested outcomes | A stated outcome for work that has not shipped. The repos' own wording is "explore whether software can make research easier to examine" |
| "the bleeding edge of Artificial Intelligence" | `lab-voice`: "Avoid ... exaggerated AI capability"; plain language | Decorative, and it raises the bar the lab then has to clear on camera |
| No "research, not medical advice" statement | `cf-research-context`: "Community-facing materials must carry the 'research, not medical advice' boundary" | Every repo README carries it. The channel is the most community-facing surface the lab has |
| No independence statement, while three of four links point at the Foundation | `cf-research-context`: "Do not ... imply affiliation." `lab-voice`: donation links "free of ... implied affiliation" | Without the disclaimer, the link block reads as an official chapter channel |

What is already right and should stay: the honest line that stream revenue funds API tokens, the clear separation of that from the CFF donation link, and the plain admission that the channel is also for WebGL, sim racing, game design, and guitar. That mix is a feature; it tells people who you are. It only needs the CF material to be clearly labelled so a clip is never mistaken for guidance.

## Proposed About text

Written to the `lab-voice` register for a community welcome. Replace the first paragraph and add the boundary block; keep the rest.

> Welcome to the Sensei Ewok Research Lab, an independent, donated-time effort. We build small open-source tools that try to make cystic fibrosis research easier to examine: tracing published numbers back to their source, checking claims against primary literature, and writing down what we learn so a stranger can verify it. Everything we make is built in public and open for inspection, including the mistakes.
>
> We are not clinicians and we do not give medical advice. Nothing here diagnoses, treats, or recommends a therapy. We are not affiliated with the Cystic Fibrosis Foundation; we link to them because they do the groundwork, and if you choose to give, that is where it should go.
>
> This channel is also a workshop for the other things I love: 3D WebGL development, sim racing, game design, and guitar when it is time to decompress. CF research streams are labelled as such.
>
> Stream revenue goes toward the API costs of keeping the lab running. To support the wider mission directly, donate to the Cystic Fibrosis Foundation (Washington Chapter): https://give.cff.org/washington/donate

Word count is close to the original. "Independent" and "not affiliated" are each said once, which is enough; repeating them reads as anxiety.

## Standing boundary for video descriptions and end cards

One block, pasted unchanged under every CF-related video. Consistency is what makes it credible.

> Research, not medical advice. The Sensei Ewok Research Lab is an independent, open-source project with no clinical role and no affiliation with the Cystic Fibrosis Foundation. Nothing in this video diagnoses, treats, or recommends a therapy. Decisions about care belong with you and your care team. Sources for every number shown: [link to the repo file for this video].

The last sentence is the one that costs effort and earns the most. It commits every video to the same evidence ledger the repos already require, and it makes the channel the public face of the catalog rather than a separate voice that has to be reconciled with it later.

## Not done

- No channel settings were changed; I have read access only, which is correct.
- The senseiewok.ai site text was not reviewed. It likely has the same first paragraph and should be read against the same table.
