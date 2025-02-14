from PIL import Image
import numpy as np
from skimage import measure, morphology
import matplotlib.pyplot as plt
import os
base_dir = os.getcwd()  # Current working directory
file_name = "Cpop_grids 2_123_17470_Cpop_grids 2.txt_monomers_placement.png"

file_path = os.path.join(base_dir, "plots", "tmp", file_name)

# Load the image and convert to grayscale
image_path = file_path  # Replace with your file path
img = Image.open(image_path).convert('L')  # Convert to grayscale
img_array = np.array(img)

# Thresholding: Identify black pixels as monomer regions
binary_image = img_array < 50  # Adjust threshold if needed

# Apply morphological dilation to connect nearby regions
dilated_image = morphology.dilation(binary_image, morphology.square(3))

# Calculate density
def calculate_density(binary_image):
    edge_pixels = np.sum(binary_image)
    total_pixels = binary_image.size
    return edge_pixels / total_pixels

# Count number of connected components (islands)
def count_islands(binary_image):
    labeled_array, num_features = measure.label(binary_image, background=0, return_num=True)
    return num_features

# Calculate the length of the longest connected component
def longest_strand_length(binary_image):
    labeled_array, num_features = measure.label(binary_image, background=0, return_num=True)
    regions = measure.regionprops(labeled_array)
    if regions:
        return max(region.area for region in regions)
    return 0

# Calculate average strand size
def average_strand_length(binary_image):
    labeled_array, num_features = measure.label(binary_image, background=0, return_num=True)
    regions = measure.regionprops(labeled_array)
    if regions:
        return np.mean([region.area for region in regions])
    return 0


# Calculate strand lengths
def get_strand_lengths(binary_image):
    labeled_array, num_features = measure.label(binary_image, background=0, return_num=True)
    regions = measure.regionprops(labeled_array)
    strand_lengths = [region.area for region in regions]
    return strand_lengths

# Get the lengths of all strands
strand_lengths = get_strand_lengths(dilated_image)

# Plot the histogram of strand lengths
plt.figure(figsize=(10, 6))
plt.hist(strand_lengths, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of Strand Lengths")
plt.xlabel("Strand Length (pixels)")
plt.ylabel("Frequency")
plt.show()

# Print summary statistics
if strand_lengths:
    print(f"Number of Strands: {len(strand_lengths)}")
    print(f"Mean Strand Length: {np.mean(strand_lengths):.2f} pixels")
    print(f"Max Strand Length: {max(strand_lengths)} pixels")
else:
    print("No strands detected.")


# Perform analysis
density = calculate_density(dilated_image)
num_islands = count_islands(dilated_image)
longest_strand = longest_strand_length(dilated_image)
average_strand = average_strand_length(dilated_image)

# Display results
print("Density:", density)
print("Number of Islands:", num_islands)
print("Longest Strand Length:", longest_strand)
print("Average Strand Size:", average_strand)

# Visualize the original and processed images
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
ax[0].imshow(binary_image, cmap='gray')
ax[0].set_title("Original Binary Image")
ax[0].axis('off')

ax[1].imshow(dilated_image, cmap='gray')
ax[1].set_title("Dilated Image")
ax[1].axis('off')

plt.tight_layout()
plt.show()
