import numpy as np
import matplotlib.pyplot as plt
import random

# Function to generate a synthetic monomer strand plot
def generate_strand_image(image_size=(600, 600), num_monomers=3000, connectivity="random", save_path=None):
    """
    Generates a small synthetic graph-like strand pattern and saves or displays the image.

    Parameters:
    - image_size: Tuple (width, height) of the image
    - num_monomers: Number of monomer points
    - connectivity: Type of connectivity ("random", "linear", "branching")
    - save_path: File path to save the generated image
    """
    img = np.ones(image_size) * 255  # White background
    x_center, y_center = image_size[0] // 2, image_size[1] // 2  # Start at center

    # Generate coordinates based on the selected connectivity type
    points = [(x_center, y_center)]
    for _ in range(num_monomers - 1):
        x, y = points[-1]

        if connectivity == "random":
            x += random.choice([-1, 0, 1]) * 5
            y += random.choice([-1, 0, 1]) * 5
        elif connectivity == "linear":
            x += 5  # Move to the right
        elif connectivity == "branching":
            if random.random() > 0.5:
                x += 5
            else:
                y += 5

        # Keep within bounds
        x = max(0, min(image_size[0] - 1, x))
        y = max(0, min(image_size[1] - 1, y))
        points.append((x, y))

    # Draw points on the image
    for x, y in points:
        img[int(y), int(x)] = 0  # Black pixel for monomers

    # Save or display the image
    if save_path:
        plt.imsave(save_path, img, cmap="gray")
        print(f"Image saved to {save_path}")
    else:
        plt.imshow(img, cmap="gray")
        plt.title(f"Synthetic Monomer Strand ({connectivity} connectivity)")
        plt.axis("off")
        plt.show()

# Generate test images
generate_strand_image(connectivity="random")
generate_strand_image(connectivity="linear")
generate_strand_image(connectivity="branching")

