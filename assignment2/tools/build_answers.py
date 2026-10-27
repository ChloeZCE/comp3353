"""Fill @@path@@ / @@path:first-last@@ markers in ANSWERS.template.md with file contents -> ANSWERS.md."""
import re
from pathlib import Path

root = Path(__file__).resolve().parent.parent
tpl = (root / "tools" / "ANSWERS.template.md").read_text()


def include(m):
    path, _, span = m.group(1).partition(":")
    text = (root / path).read_text().rstrip("\n")
    if span:
        a, b = map(int, span.split("-"))
        text = "\n".join(text.split("\n")[a - 1:b])
    return text


out = re.sub(r"@@(.+?)@@", include, tpl)
(root / "ANSWERS.md").write_text(out)

# Word/PDF versions for Moodle:
#   pandoc ANSWERS.md -o ANSWERS.docx --resource-path=.
#   soffice --headless --convert-to pdf ANSWERS.docx
