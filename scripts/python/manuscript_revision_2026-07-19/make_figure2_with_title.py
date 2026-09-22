import fitz  # PyMuPDF
from PIL import Image, ImageDraw, ImageFont
import io, os, sys

# Paths
fig_png_path = r'E:/workbuddy/Genetic Epidemiology/投稿资料/03_主文图_Figures_Main/Figure2.png'
title_img_path = r'C:/Users/Administrator/.workbuddy/clipboard-images/clipboard-2026-07-19T11-40-25-246Z-f9da8002.png'
out_pdf_path = r'E:/workbuddy/Genetic Epidemiology/投稿资料/03_主文图_Figures_Main/Figure2.pdf'

# 1. Load Figure2.png
fig = Image.open(fig_png_path).convert('RGB')
fig_w, fig_h = fig.size
print('Figure2.png size:', fig_w, fig_h)

# 2. Detect and crop the top title from Figure2.png
# Find horizontal projection (non-white pixels per row)
gray = fig.convert('L')
px = gray.load()
proj_y = [0] * fig_h
for y in range(fig_h):
    s = 0
    for x in range(fig_w):
        s += 255 - px[x, y]
    proj_y[y] = s

# Normalize and find the bottom of the title (a long stretch of low projection after title)
max_p = max(proj_y)
# Find first major drop below threshold, followed by sustained low region (title-content gap)
# The title is at the top, then a gap, then the plot body.
# Find the bottom of the title by locating the longest low-projection stretch at the top.
th = max_p * 0.02
best_start = 0
best_len = 0
best_end = 0
cur_start = 0
cur_len = 0
for y in range(fig_h):
    if proj_y[y] < th:
        if cur_len == 0:
            cur_start = y
        cur_len += 1
    else:
        if cur_len > best_len:
            best_len = cur_len
            best_start = cur_start
            best_end = y
        cur_len = 0
if cur_len > best_len:
    best_len = cur_len
    best_start = cur_start
    best_end = fig_h

if best_len > 20 and best_end < 300:
    crop_y = best_end
    print('Long blank gap detected:', best_start, '-', best_end, '=> crop_y =', crop_y)
else:
    crop_y = 150
    print('No clear blank gap; using crop_y =', crop_y)

cropped = fig.crop((0, crop_y, fig_w, fig_h))
print('Cropped plot size:', cropped.size)

# 3. Load title image and modify "Figure 1." -> "Figure 2."
title = Image.open(title_img_path).convert('RGB')
# Find '1' position by projections
gray_t = title.convert('L')
px_t = gray_t.load()
w, h = title.size
# vertical projection to find segments
proj = [0] * w
for x in range(w):
    s = 0
    for y in range(h):
        s += 255 - px_t[x, y]
    proj[x] = s
maxp = max(proj)
th = maxp * 0.05
segments = []
in_seg = False
start = 0
for x in range(w):
    if proj[x] > th and not in_seg:
        start = x
        in_seg = True
    elif proj[x] <= th and in_seg:
        segments.append((start, x))
        in_seg = False
if in_seg:
    segments.append((start, w))

# horizontal projection for y range
proj_y = [0] * h
for y in range(h):
    s = 0
    for x in range(w):
        s += 255 - px_t[x, y]
    proj_y[y] = s
yth = max(proj_y) * 0.05
y_start = next(y for y in range(h) if proj_y[y] > yth)
y_end = next(y for y in range(h-1, -1, -1) if proj_y[y] > yth)
char_h = y_end - y_start + 1

print('Title segments:', segments[:10])
print('y_start', y_start, 'y_end', y_end, 'char_h', char_h)

if len(segments) < 8:
    print('ERROR: could not locate "1" in title image', segments)
    sys.exit(1)

one_x1, one_x2 = segments[6]
dot_x1, dot_x2 = segments[7]
print('Detected "1" segment:', one_x1, one_x2, 'dot:', dot_x1, dot_x2)

# Cover the '1' with white
draw = ImageDraw.Draw(title)
erase_x1 = one_x1 - 2
erase_x2 = dot_x1 - 1
draw.rectangle([erase_x1, y_start, erase_x2, y_end], fill=(255, 255, 255))

# Draw '2' matching size
font_candidates = [
    r'C:/Windows/Fonts/times.ttf',
    r'C:/Windows/Fonts/arial.ttf',
    r'C:/Windows/Fonts/calibri.ttf',
]
font_path = None
for p in font_candidates:
    if os.path.exists(p):
        font_path = p
        break

text = '2'
best_font = None
best_size = 0
best_err = 1e9
if font_path:
    for size in range(1, 80):
        try:
            fnt = ImageFont.truetype(font_path, size)
            bbox = draw.textbbox((0, 0), text, font=fnt)
            h2 = bbox[3] - bbox[1]
            err = abs(h2 - char_h)
            if err < best_err:
                best_err = err
                best_size = size
                best_font = fnt
        except Exception:
            continue
if best_font is None:
    best_font = ImageFont.load_default()
    bbox = draw.textbbox((0, 0), text, font=best_font)
else:
    bbox = draw.textbbox((0, 0), text, font=best_font)

text_w = bbox[2] - bbox[0]
text_h = bbox[3] - bbox[1]
print('Best font size for "2":', best_size, 'text_h', text_h, 'char_h', char_h)

x = erase_x1 + (erase_x2 - erase_x1 - text_w) // 2
y = (y_start + y_end - text_h) // 2 - bbox[1]
draw.text((x, y), text, fill=(0, 0, 0), font=best_font)
print('Drew "2" at', (x, y))

title.save('title_modified_preview.png')
print('saved modified title preview to title_modified_preview.png')

# 4. Resize title image to match figure width
title_w = fig_w
title_h = int(title.size[1] * fig_w / title.size[0])
title_resized = title.resize((title_w, title_h), Image.LANCZOS)

# 5. Concatenate vertically (plot on top, title below)
combined = Image.new('RGB', (fig_w, cropped.size[1] + title_h), (255, 255, 255))
combined.paste(cropped, (0, 0))
combined.paste(title_resized, (0, cropped.size[1]))

# 6. Save as PDF and PNG
combined.save(out_pdf_path, 'PDF', resolution=200.0)
print('Saved new Figure2.pdf to', out_pdf_path)
combined.save(r'E:/workbuddy/Genetic Epidemiology/投稿资料/03_主文图_Figures_Main/Figure2.png', 'PNG')
print('Saved new Figure2.png to same directory')
combined.save('Figure2_combined_preview.png', 'PNG')
print('Saved combined preview to Figure2_combined_preview.png')
