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
    total_width = sum(widths)
    max_height = max(heights)

    # Create an empty image with the combined size
    panel = Image.new("RGBA", (total_width, max_height))

    # Paste images onto the panel
    x_offset = 0
    for img in images:
        panel.paste(img, (x_offset, 0))
        x_offset += img.size[0]

    # Save the panel image in the subfolder
    panel.save(os.path.join(subfolder, "panel.png"))

if __name__ == "__main__":
    subfolder = sys.argv[1]
    create_panel(subfolder)
