"""Fill @@path@@ markers in ANSWERS.template.md with file contents -> ANSWERS.md."""
import re
from pathlib import Path

root = Path(__file__).resolve().parent.parent
tpl = (root / "tools" / "ANSWERS.template.md").read_text()
out = re.sub(r"@@(.+?)@@", lambda m: (root / m.group(1)).read_text().rstrip("\n"), tpl)
(root / "ANSWERS.md").write_text(out)
