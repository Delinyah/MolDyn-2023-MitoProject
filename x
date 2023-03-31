    font = ImageFont.truetype("arial.ttf", 24)
    text_width, text_height = draw.textsize(title_text, font=font)
    draw.rectangle([(0, panel.height - text_height - 10), (text_width + 20, panel.height)], fill=(200, 200, 200), outline=None)
    draw.text((10, panel.height - 30), title_text, font=font, fill=(255, 0, 0))
