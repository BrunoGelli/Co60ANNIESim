"""Read flat RATDB blocks for static checks; deliberately not a RAT resolver.

Reject unsupported syntax. Runtime DB precedence still requires the tuned RAT core.
"""
import json
import re


def blocks(path):
    text = path.read_text()
    # Remove comments outside quoted strings only.
    text = re.sub(r'"(?:\\.|[^"\\])*"|//[^\n]*|/\*[\s\S]*?\*/',
                  lambda m: m[0] if m[0].startswith('"') else '', text)
    result = []
    for body in re.findall(r'\{([^{}]*)\}', text):
        body = re.sub(r'\b([A-Za-z_][A-Za-z_0-9]*)\s*:', r'"\1":', body)
        # RAT accepts e.g. 00.0; JSON does not. Preserve quoted strings exactly.
        body = re.sub(r'"(?:\\.|[^"\\])*"|(?<![\w.])[+-]?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?',
                      lambda m: m[0] if m[0].startswith('"') else
                      json.dumps(float(m[0]) if any(c in m[0] for c in '.eE') else int(m[0])), body)
        body = re.sub(r',\s*$', '', body)
        body = re.sub(r',\s*\]', ']', body)
        result.append(json.loads('{' + body + '}'))
    return result


def dump_blocks(records):
    return '\n'.join(json.dumps(b, indent=2) for b in records) + '\n'
