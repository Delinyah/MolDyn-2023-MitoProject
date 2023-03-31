import os
import sys
from PIL import Image, ImageDraw, ImageFont

def create_panel(subfolder):
    # List of image file names in the order you want them to appear in the panel
    image_files = [
        "Membrane_Protein_Angles.png",
        "Log1-1_enrichment.png",
        "Log1-1_enrichment_hm.png",
        "Cormat_lip_comp.png",
    ]

    # Read images from the subfolder
    images = [Image.open(os.path.join(subfolder, file)) for file in image_files]

    # Set the target height for all images
    target_height = min(img.height for img in images)

    # Resize images to the same height while maintaining the aspect ratio
    for i, img in enumerate(images):
        img.thumbnail((img.width, target_height), Image.ANTIALIAS)
        images[i] = img

    # Calculate the size of the combined image
    widths, heights = zip(*(img.size for img in images))
    total_width = sum(widths)
    max_height = max(heights)

    # Create an empty image with the combined size
    title_width = 200
    panel = Image.new("RGBA", (total_width + title_width, max_height))

    # Draw the folder name as a title on the left side of the panel
    draw = ImageDraw.Draw(panel)
    title_text = os.path.basename(subfolder)
    font = ImageFont.truetype("arial.ttf", 24)
    draw.text((10, 10), title_text, font=font, fill=(0, 0, 0))

    # Paste the title image onto the top left corner of the panel
    title_image = Image.new("RGBA", (title_width, max_height), color=(255, 255, 255))
    title_draw = ImageDraw.Draw(title_image)
    title_draw.text((10, 10), "Panel", font=font, fill=(0, 0, 0))
    panel.paste(title_image, (0, 0))

    # Paste images onto the panel
    x_offset = title_width
    space_between_images = 50 # number of pixels to move images closer
    for img in images:
        panel.paste(img, (x_offset, 0))
        x_offset += img.size[0] - space_between_images

    # Save the panel image in the subfolder
    panel.save(os.path.join(subfolder, "panel.png"))

if __name__ == "__main__":
    subfolder = sys.argv[1]
    create_panel(subfolder)
