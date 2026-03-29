# Load the dplyr package for data manipulation functions like filter(), mutate(), and arrange()
library(dplyr)

# Load the ggplot2 package for creating plots
library(ggplot2)

# Load the ggrepel package to add non-overlapping text labels to plots
library(ggrepel)

# Set the working directory to the folder where the input file is located
setwd("/Users/sindypinero/Downloads/")

# Read the CSV file into an object called df
# header = TRUE means the first row contains column names
# stringsAsFactors = FALSE keeps text columns as character strings instead of factors
df <- read.csv("DCE_results_1.csv", header = TRUE, stringsAsFactors = FALSE)

# Keep only rows where predictor and outcome are not equal to "#N/A"
df <- df %>%
  filter(predictor != "#N/A", outcome != "#N/A")

# Create a new data frame called df_top
df_top <- df %>%
  mutate(
    # Create a new column with the absolute value of DCE
    abs_dce    = abs(DCE),
    
    # Create a label combining predictor and outcome with "->" between them
    path_label = paste(predictor, outcome, sep = "->")
  ) %>%
  # Sort the rows from the largest to the smallest absolute DCE
  arrange(desc(abs_dce))

# Show summary statistics for the abs_dce column
summary(df_top$abs_dce)

# Show the 75th, 90th, and 95th percentiles of the abs_dce column
quantile(df_top$abs_dce, probs = c(0.75, 0.90, 0.95))

# Create a dot plot object and store it in dot_plot
dot_plot <- ggplot(df_top, aes(x = p_value, y = abs_dce)) +
  
  # Add points to the plot
  # x-axis = p_value, y-axis = abs_dce, and point color represents DCE
  geom_point(aes(color = DCE), size = 2) +
  
  # Add text labels to selected points, using ggrepel to avoid label overlap
  geom_text_repel(
    # Only label rows where adjusted p-value is < 0.05 and abs_dce is > 0.41
    data = subset(df_top, adj_p_value < 0.05 & abs_dce > 0.41),
    
    # Use the path_label column as the text label
    aes(label = path_label),
    
    # Set the text size
    size = 3,
    
    # Allow all selected labels to be plotted, even if there are many
    max.overlaps = Inf,
    
    # Increase repulsion force between labels
    force = 5,
    
    # Add padding around the label boxes
    box.padding = 0.4,
    
    # Add padding around the data points
    point.padding = 0.3,
    
    # Set line thickness for the connecting segments
    segment.size = 0.3,
    
    # Set segment color to grey
    segment.color = "grey50",
    
    # Set a seed so label placement is reproducible
    seed = 42
  ) +
  
  # Use a diverging color scale:
  # blue for negative DCE, white around zero, red for positive DCE
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "DCE") +
  
  # Show the x-axis on a log10 scale
  scale_x_continuous(trans = "log10") +
  
  # Set axis labels
  labs(x = "p-value (log scale)", y = "Absolute DCE") +
  
  # Apply a minimal theme with base font size 11
  theme_minimal(base_size = 11) +
  
  # Adjust the size of y-axis text labels
  theme(axis.text.y = element_text(size = 10))

# Display the plot in the R plotting window
print(dot_plot)

# Save the plot to a PNG file
ggsave(
  # Name of the output file
  filename = "Gene_pairs_DCE_dotplot.png",
  
  # Plot object to save
  plot     = dot_plot,
  
  # File format/device
  device   = "png",
  
  # Width of the saved image in inches
  width    = 10,
  
  # Height of the saved image in inches
  height   = 8,
  
  # Resolution in dots per inch
  dpi      = 300
)

# Print a message showing the folder where the file was saved
message("Saved to: ", getwd())