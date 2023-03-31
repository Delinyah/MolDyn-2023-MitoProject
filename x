import os
import sys
from PIL import Image

def create_panel(subfolder):
    # List of image file names in the order you want them to appear in the panel
    image_files = [
        "Membrane_Protein_Angles.png",
        "Log1-1_enrichment.png",
        "Log1-1_enrichment_hm.png",
        "Covmat_lip_comp,.png",
    ]

    # Read images from the subfolder
    images = [Image.open(os.path.join(subfolder, file)) for file in image_files]

    # Calculate the size of the combined image
    widths, heights = zip(*(img.size for img in images))
    total_height = sum(heights)
    max_width = max(widths)

    # Create an empty image with the combined size
    panel = Image.new("RGBA", (max_width, total_height))

    # Paste images onto the panel
    y_offset = 0
    for img in images:
        panel.paste(img, (0, y_offset))
        y_offset += img.size[1]

    # Save the panel image in the subfolder
    panel.save(os.path.join(subfolder, "panel.png"))

if __name__ == "__main__":
    subfolder = sys.argv[1]
    create_panel(subfolder)
