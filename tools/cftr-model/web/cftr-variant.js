// Reads one protein-level variant name and returns the residue numbers it names. Pure: no page, no data, no network. The page (cftr-page.js) checks the
// reference letters against the UniProt sequence and says nothing about what a change does.
// Accepted: F508del, delF508, G551D, W1282X / W1282* / W1282Ter, I507_F508del, the same with three-letter codes and an optional "p." and brackets
// (p.Phe508del, p.(Gly551Asp), p.Ile507_Phe508del), and a bare residue number. Refused, with a reason: DNA-level and intron names, frameshifts, more than one name.
const AA = 'ACDEFGHIKLMNPQRSTVWY';
const AA3 = { Ala: 'A', Arg: 'R', Asn: 'N', Asp: 'D', Cys: 'C', Gln: 'Q', Glu: 'E', Gly: 'G', His: 'H', Ile: 'I', Leu: 'L', Lys: 'K', Met: 'M', Phe: 'F', Pro: 'P', Ser: 'S', Thr: 'T', Trp: 'W', Tyr: 'Y', Val: 'V' };
const R = '([a-z]{3}|[a-z])(\\d{1,5})';
const FORMS = ['^del' + R + '$', '^' + R + '_' + R + 'del$', '^' + R + 'del$', '^' + R + '([a-z]{3}|[a-z]|\\*)$'].map(f => new RegExp(f, 'i'));

const letter = (t) => (t.length === 3 ? AA3[t[0].toUpperCase() + t.slice(1).toLowerCase()] : AA.includes(t.toUpperCase()) ? t.toUpperCase() : undefined);

export function parseVariant(text, max = 1480) {
  const t = String(text ?? '').trim(), no = (why) => ({ ok: false, why });
  if (!t) return no('empty');
  if (/^\d+$/.test(t)) { const n = Number(t); return n >= 1 && n <= max ? { ok: true, start: n, end: n, ref: [], kind: 'number' } : no('range'); }
  if (/^[cgnrm]\.|^ivs|\d[+-]\d|>|^\d+(del|ins|dup)/i.test(t)) return no('dna');
  if (/[\s/;,|&+]/.test(t)) return no('many');
  if (/fs/i.test(t)) return no('fs');
  const s = t.replace(/^p\./i, '').replace(/^\((.*)\)$/, '$1');
  const [i, m] = FORMS.map((re, k) => [k, s.match(re)]).find(([, x]) => x) || [];
  if (!m) return no('unread');
  const range = i === 1, start = Number(m[2]), end = range ? Number(m[4]) : start;
  const ref = range ? [letter(m[1]), letter(m[3])] : [letter(m[1])];
  const to = i === 3 ? (/^(x|\*|ter)$/i.test(m[3]) ? '*' : letter(m[3])) : '-';
  if (ref.includes(undefined) || !to) return no('unread');
  if (start < 1 || end > max || (range && end <= start)) return no('range');
  return { ok: true, start, end, ref, kind: i < 3 ? 'del' : to === '*' ? 'stop' : 'sub' };
}
