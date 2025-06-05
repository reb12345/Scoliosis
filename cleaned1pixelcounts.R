# ==============================================================================
# SPINE IMAGE ANALYSIS - DICOM POST-PROCESSING
# Analysis of cumulative curvature measurements from Frontera 820 and 960 datasets
# ==============================================================================

# Required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# ==============================================================================
# DATA INITIALIZATION AND PREPROCESSING
# ==============================================================================

# Initialize datasets from method 8 results
# Trial 6 cumulative type processing for both 820 and 960 resolution datasets
frontera820_method2_trial6_interpret1 <- fr0ntera820_1method8
frontera820_method2_trial6_interpret2 <- fr0ntera820_2method8
frontera960_method2_trial6_interpret1 <- fr0ntera960_1method8
frontera960_method2_trial6_interpret2 <- fr0ntera960_2method8

# ==============================================================================
# DATA CLEANING - REMOVE HEADER ROWS
# ==============================================================================

# Remove first 2 rows (likely headers) from each dataset
# Note: Different end points suggest varying data sizes
frontera820_method2_trial6_interpret1 <- slice(frontera820_method2_trial6_interpret1, 3:10991)
frontera820_method2_trial6_interpret2 <- slice(frontera820_method2_trial6_interpret2, 3:10992)

frontera960_method2_trial6_interpret1 <- slice(frontera960_method2_trial6_interpret1, 3:8736)
frontera960_method2_trial6_interpret2 <- slice(frontera960_method2_trial6_interpret2, 3:8731)

# ==============================================================================
# DATASET COMBINATION BY RESOLUTION
# ==============================================================================

# Combine interpretation datasets by resolution type
frontera820_method2_interpretations <- rbind(frontera820_method2_trial6_interpret1, 
                                             frontera820_method2_trial6_interpret2)

frontera960_method2_interpretations <- rbind(frontera960_method2_trial6_interpret1, 
                                             frontera960_method2_trial6_interpret2)

# ==============================================================================
# DATA SELECTION AND MASTER DATASET CREATION
# ==============================================================================

# Select relevant columns for analysis (image_name through cumulative2)
frontera820_method2_interpretations_join <- frontera820_method2_interpretations %>% 
  select(image_name:cumulative2)

frontera960_method2_interpretations_join <- frontera960_method2_interpretations %>% 
  select(image_name:cumulative2)

# Create initial combined dataset
frontera_method2_interpretations <- rbind(frontera960_method2_interpretations_join,
                                          frontera820_method2_interpretations_join)

# Add resolution identifier flags for downstream analysis
frontera820_method2_interpretations <- frontera820_method2_interpretations %>% 
  mutate(is820 = TRUE)

frontera960_method2_interpretations <- frontera960_method2_interpretations %>% 
  mutate(is820 = FALSE)

# Create final combined dataset with resolution flags
frontera_method2_interpretations <- rbind(frontera820_method2_interpretations,
                                          frontera960_method2_interpretations)

# ==============================================================================
# FRONTERA 960 PROCESSING AND CALIBRATION
# ==============================================================================

# Apply calibration factors to 960 dataset measurements
frontera960_method2_interpretations <- frontera960_method2_interpretations %>% 
  left_join(img_beta_df_1_) %>%
  # Apply calibration Value to all measurement variables
  mutate(
    curvature_method_2 = curvature_method_2 * Value,
    cumulative1 = cumulative1 * Value,
    cumulative2 = cumulative2 * Value,
    curvature_method_3 = curvature_method_3 * Value
  ) %>%
  # Remove records without calibration values
  filter(!is.na(Value))

# Add patient EID information
frontera960_method2_interpretations <- left_join(frontera960_method2_interpretations, 
                                                 Patient_EID)

# Apply calibration to X960 dataset (legacy processing)
X960_method2_interpretations <- X960_method2_interpretations %>% 
  mutate(
    cumulative1 = cumulative1 * Value,
    cumulative2 = cumulative2 * Value
  )

# ==============================================================================
# DATA QUALITY CONTROL - REMOVE PROBLEMATIC IMAGES
# ==============================================================================

# Define problematic image names for 960 dataset
problematic_images_960 <- c(
  "1.2.840.113619.2.110.210419.20171120183234.1.12.12.1",
  "1.2.840.113619.2.110.210419.20171128154331.1.12.12.1"
)

# Remove problematic images from 960 dataset
frontera960_method2_interpretations <- frontera960_method2_interpretations %>% 
  filter(!str_detect(image_name, paste(problematic_images_960, collapse = "|")))

# ==============================================================================
# DUPLICATE ANALYSIS - FRONTERA 960
# ==============================================================================

# Identify patients with multiple measurements
duplicate_eids_960 <- names(table(frontera960_method2_interpretations$eid))[
  table(frontera960_method2_interpretations$eid) > 1
]

# Create comparison dataset for duplicates
duplicate_comparison_960 <- frontera960_method2_interpretations[
  frontera960_method2_interpretations$eid %in% duplicate_eids_960, 
] %>%
  group_by(eid) %>% 
  summarise(
    measurement_combined = list(cumulative2),
    image_name_combined = list(image_name)
  ) %>%
  # Parse combined measurements into separate columns
  separate(measurement_combined, sep = ",", into = c("first", "second")) %>%
  separate(image_name_combined, sep = ",", into = c("first_image", "second_image"))

# Clean and convert measurement values
duplicate_comparison_960 <- duplicate_comparison_960 %>%
  mutate(
    first = as.numeric(str_remove_all(first, "c\\(")),
    second = as.numeric(str_remove_all(second, "\\)")),
    difference = abs(first - second)
  ) %>%
  na.omit()

# Assess correlation between duplicate measurements
correlation_960_duplicates <- cor(duplicate_comparison_960$first, 
                                  duplicate_comparison_960$second)
print(paste("Correlation between duplicate 960 measurements:", 
            round(correlation_960_duplicates, 3)))

# Visualize duplicate measurement agreement
duplicate_plot_960 <- ggplot(duplicate_comparison_960, 
                             aes(x = first, y = second, color = difference)) +
  geom_point() +
  labs(
    title = "Duplicate Measurement Agreement - Frontera 960",
    x = "First Measurement (cumulative2)",
    y = "Second Measurement (cumulative2)",
    color = "Absolute Difference"
  ) +
  theme_minimal()

print(duplicate_plot_960)

# ==============================================================================
# FRONTERA 820 PROCESSING AND CALIBRATION
# ==============================================================================

# Apply calibration factors to 820 dataset measurements
frontera820_method2_interpretations <- frontera820_method2_interpretations %>% 
  left_join(img_beta_df_1_) %>%
  # Apply calibration Value to all measurement variables
  mutate(
    curvature_method_2 = curvature_method_2 * Value,
    cumulative1 = cumulative1 * Value,
    cumulative2 = cumulative2 * Value,
    curvature_method_3 = curvature_method_3 * Value
  ) %>%
  # Remove records without calibration values
  filter(!is.na(Value))

# Add patient EID information
frontera820_method2_interpretations <- left_join(frontera820_method2_interpretations, 
                                                 Patient_EID)

# ==============================================================================
# DATA QUALITY CONTROL - REMOVE PROBLEMATIC IMAGES (820)
# ==============================================================================

# Define problematic image names for 820 dataset
problematic_images_820 <- c(
  "1.2.840.113619.2.110.212038.20180426143443.1.12.12.1",
  "1.2.840.113619.2.110.212038.20190924115521.1.3.12.1"
)

# Remove problematic images from 820 dataset
frontera820_method2_interpretations <- frontera820_method2_interpretations %>% 
  filter(!str_detect(image_name, paste(problematic_images_820, collapse = "|")))

# ==============================================================================
# DUPLICATE ANALYSIS - FRONTERA 820
# ==============================================================================

# Identify patients with multiple measurements
duplicate_eids_820 <- names(table(frontera820_method2_interpretations$eid))[
  table(frontera820_method2_interpretations$eid) > 1
]

# Create comparison dataset for duplicates (using cumulative1 for 820)
duplicate_comparison_820 <- frontera820_method2_interpretations[
  frontera820_method2_interpretations$eid %in% duplicate_eids_820, 
] %>%
  group_by(eid) %>% 
  summarise(
    measurement_combined = list(cumulative1),
    image_name_combined = list(image_name)
  ) %>%
  # Parse combined measurements into separate columns
  separate(measurement_combined, sep = ",", into = c("first", "second")) %>%
  separate(image_name_combined, sep = ",", into = c("first_image", "second_image"))

# Clean and convert measurement values
duplicate_comparison_820 <- duplicate_comparison_820 %>%
  mutate(
    first = as.numeric(str_remove_all(first, "c\\(")),
    second = as.numeric(str_remove_all(second, "\\)")),
    difference = abs(first - second)
  ) %>%
  na.omit()

# Assess correlation between duplicate measurements
correlation_820_duplicates <- cor(duplicate_comparison_820$first, 
                                  duplicate_comparison_820$second)
print(paste("Correlation between duplicate 820 measurements:", 
            round(correlation_820_duplicates, 3)))

# Visualize duplicate measurement agreement
duplicate_plot_820 <- ggplot(duplicate_comparison_820, 
                             aes(x = first, y = second, color = difference)) +
  geom_point() +
  labs(
    title = "Duplicate Measurement Agreement - Frontera 820",
    x = "First Measurement (cumulative1)",
    y = "Second Measurement (cumulative1)",
    color = "Absolute Difference"
  ) +
  theme_minimal()

print(duplicate_plot_820)

# ==============================================================================
# IMAGE ANALYSIS INTEGRATION
# ==============================================================================

# Clean image names for joining (remove .png extension)
spine_image_analysis$image_name <- str_remove_all(spine_image_analysis$image_name, ".png")

# Add image analysis metrics to combined dataset
frontera_method2_interpretations <- frontera_method2_interpretations %>%
  left_join(spine_image_analysis) %>%
  # Calculate white to black pixel ratio
  mutate(white_black_ratio = num_white_pixels / num_black_pixels)

# Generate exploratory histograms for image analysis metrics
histogram_vars <- c("num_white_pixels", "num_black_pixels", "min_x", "min_y", 
                    "max_x", "max_y", "white_black_ratio")

for (var in histogram_vars) {
  hist_plot <- hist(frontera_method2_interpretations[[var]], 
                    breaks = 50, 
                    main = paste("Distribution of", var),
                    xlab = var)
}

# Add image analysis to individual resolution datasets
spine_image_analysis <- spine_image_analysis %>% 
  mutate(white_black_ratio = num_white_pixels / num_black_pixels)

frontera960_method2_interpretations <- left_join(frontera960_method2_interpretations, 
                                                 spine_image_analysis)
frontera820_method2_interpretations <- left_join(frontera820_method2_interpretations, 
                                                 spine_image_analysis)

# ==============================================================================
# MORPHOMETRIC CALCULATIONS
# ==============================================================================

# Calculate spine height from bounding box coordinates
frontera_method2_interpretations <- frontera_method2_interpretations %>% 
  mutate(
    calculated_spine_height = sqrt((min_x - max_x)^2 + (min_y - max_y)^2)
  ) %>%
  # Join calibration data and apply to spine height
  left_join(fid12144) %>%
  mutate(calculated_spine_height = calculated_spine_height * Value)

# Parse coordinate points and calculate rotated spine length
frontera_method2_interpretations <- frontera_method2_interpretations %>%
  # Clean coordinate formatting
  mutate(
    finalTopPoint = gsub("\\[|\\]", "", finalTopPoint),
    finalBottomPoint = gsub("\\[|\\]", "", finalBottomPoint)
  ) %>%
  # Separate coordinates into individual columns
  separate(finalTopPoint, 
           into = c("finalTopPointX", "finalTopPointY"), 
           sep = ", ", 
           convert = TRUE) %>%
  separate(finalBottomPoint, 
           into = c("finalBottomPointX", "finalBottomPointY"), 
           sep = ", ", 
           convert = TRUE) %>%
  # Calculate rotated spine length
  mutate(
    rotated_spine_length = sqrt((finalBottomPointX - finalTopPointX)^2 + 
                                  (finalBottomPointY - finalTopPointY)^2),
    rotated_spine_length = rotated_spine_length * Value
  )

# ==============================================================================
# TRIAL 8 ANALYSIS (820 DATASET SUBSET)
# ==============================================================================

# Process trial 8 data for 820 dataset
frontera820_method2_trial8_interpret1 <- slice(frontera820_method2_trial8_interpret1, 3:10991) %>%
  # Apply calibration
  left_join(img_beta_df_1_) %>%
  mutate(
    curvature_method_2 = curvature_method_2 * Value,
    cumulative1 = cumulative1 * Value,
    cumulative2 = cumulative2 * Value,
    curvature_method_3 = curvature_method_3 * Value
  ) %>%
  filter(!is.na(Value)) %>%
  # Add patient information
  left_join(Patient_EID) %>%
  # Remove problematic images
  filter(!str_detect(image_name, paste(problematic_images_820, collapse = "|"))) %>%
  # Add image analysis
  left_join(spine_image_analysis) %>%
  # Calculate morphometric measures
  mutate(
    calculated_spine_height = sqrt((min_x - max_x)^2 + (min_y - max_y)^2),
    eid = as.numeric(as.character(eid))
  ) %>%
  left_join(fid12144) %>%
  mutate(calculated_spine_height = calculated_spine_height * Value) %>%
  # Parse and calculate rotated spine measurements
  mutate(
    finalTopPoint = gsub("\\[|\\]", "", finalTopPoint),
    finalBottomPoint = gsub("\\[|\\]", "", finalBottomPoint)
  ) %>%
  separate(finalTopPoint, 
           into = c("finalTopPointX", "finalTopPointY"), 
           sep = ", ", 
           convert = TRUE) %>%
  separate(finalBottomPoint, 
           into = c("finalBottomPointX", "finalBottomPointY"), 
           sep = ", ", 
           convert = TRUE) %>%
  mutate(
    rotated_spine_length = sqrt((finalBottomPointX - finalTopPointX)^2 + 
                                  (finalBottomPointY - finalTopPointY)^2),
    rotated_spine_length = rotated_spine_length * Value,
    # Calculate original spine length
    original_spine_length = sqrt((ogbottomx - ogtopx)^2 + (ogbottomy - ogtopy)^2),
    original_spine_length = original_spine_length * Value
  )

# ==============================================================================
# FINAL ANALYSIS VARIABLES
# ==============================================================================

# Create averaged cumulative measure
frontera_method2_interpretations <- frontera_method2_interpretations %>% 
  mutate(average_cumulative = (cumulative1 + cumulative2) / 2)

# ==============================================================================
# GWAS DATASET PREPARATION
# ==============================================================================

# Create dataset for genome-wide association studies
gwas_cumulative_dataset <- frontera_method2_interpretations %>% 
  select(image_name, cumulative1, cumulative2, average_cumulative) %>%
  # Create top 10th percentile indicators for phenotype analysis
  mutate(
    top10_cumulative1 = cumulative1 >= quantile(cumulative1, 0.9, na.rm = TRUE),
    top10_cumulative2 = cumulative2 >= quantile(cumulative2, 0.9, na.rm = TRUE),
    top10_average_cumulative = average_cumulative >= quantile(average_cumulative, 0.9, na.rm = TRUE)
  )

# ==============================================================================
# DICOM METADATA PROCESSING
# ==============================================================================

# Process DICOM image metadata for path and naming consistency
dicom_img_metadata <- dicom_img_metadata %>%
  # Extract image name from full path (remove directory structure and .dcm extension)
  mutate(image_name_slimmed = str_replace(image_name, ".*/.*/(.*)\\.dcm$", "\\1")) %>%
  # Rename columns for clarity
  rename(
    image_path = image_name,
    image_name = image_name_slimmed
  ) %>%
  # Remove header row(s)
  slice(2:84806)

# Join DICOM metadata with existing interpretations for validation
frontera_method2_interpretations <- frontera_method2_interpretations %>%
  left_join(dicom_img_metadata, by = "image_name") %>%
  # Filter for records with valid pixel spacing data
  filter(!is.na(pixel_spacing_cm_X.y))

# Analyze unique image sizes in the dataset
unique_image_sizes <- frontera_method2_interpretations %>%
  distinct(origin_image_size, .keep_all = TRUE)

# ==============================================================================
# FULL BODY TRANSPARENT PATH PROCESSING
# ==============================================================================

# Process full body transparent DXA paths
full_body_transparent_path <- full_body_transparent_path %>%
  mutate(image_name_slimmed = str_replace(path_full_body_transparent_dxa, 
                                          ".*/.*/(.*)\\.dcm$", "\\1"))

# Join path information with metadata
dicom_img_metadata <- dicom_img_metadata %>%
  rename(path_full_body_transparent_dxa = image_name) %>%
  left_join(full_body_transparent_path, by = "path_full_body_transparent_dxa")

# Check unique path counts
unique_path_count <- full_body_transparent_path %>%
  summarise(unique_strings = n_distinct(path_full_body_transparent_dxa))
print(paste("Unique full body transparent paths:", unique_path_count$unique_strings))

# ==============================================================================
# DATA DEDUPLICATION AND FILTERING
# ==============================================================================

# Remove duplicate paths and filter for valid data
dicom_img_metadata <- dicom_img_metadata %>%
  distinct(path_full_body_transparent_dxa, .keep_all = TRUE)

full_body_transparent_path <- full_body_transparent_path %>%
  distinct(path_full_body_transparent_dxa, .keep_all = TRUE) %>%
  filter(!is.na(path_full_body_transparent_dxa))

# Prepare for new processing by identifying unprocessed images
dicom_img_metadata <- dicom_img_metadata %>%
  rename(image_name = image_name_slimmed.x)

# Find images not yet processed (anti-join)
unprocessed_dicom_images <- dicom_img_metadata %>%
  anti_join(frontera_method2_interpretations, by = "image_name") %>%
  # Filter for appropriate pixel spacing (< 0.5 cm)
  filter(pixel_spacing_cm_X < 0.5)

# Save current interpretations before processing new data
frontera_method2_interpretations_backup <- frontera_method2_interpretations

# ==============================================================================
# DICOM MEASUREMENTS PROCESSING
# ==============================================================================

# Process new DICOM measurements from different resolution datasets
# Extract image names from full paths
frontera820_1_dicom_measurements <- frontera820_1_dicom_measurements %>%
  mutate(image_name = str_extract(image_name, "[^/]+$")) %>%
  slice(3:17403) %>%
  mutate(is820 = TRUE)

frontera820_2_dicom_measurements <- frontera820_2_dicom_measurements %>%
  slice(3:17402) %>%
  mutate(is820 = TRUE)

frontera960_dicom_measurements <- frontera960_dicom_measurements %>%
  slice(3:1026) %>%
  mutate(is820 = FALSE) %>%
  # Add prefix to distinguish from 820 images
  mutate(image_name = str_c("1", image_name))

# Combine all DICOM measurements
dicom_measurements_combined <- bind_rows(
  frontera820_1_dicom_measurements,
  frontera820_2_dicom_measurements,
  frontera960_dicom_measurements
) %>%
  left_join(dicom_img_metadata, by = "image_name") %>%
  distinct()

# ==============================================================================
# MASTER DATASET CREATION
# ==============================================================================

# Combine all measurements (original interpretations + new DICOM measurements)
all_measurements <- bind_rows(
  frontera_method2_interpretations_backup,
  dicom_measurements_combined
) %>%
  left_join(dicom_img_metadata, by = "image_name") %>%
  distinct() %>%
  # Apply pixel spacing calibration to all measurement variables
  mutate(
    curvature_method_2 = curvature_method_2 * pixel_spacing_cm_X,
    cumulative1 = cumulative1 * pixel_spacing_cm_X,
    cumulative2 = cumulative2 * pixel_spacing_cm_X,
    curvature_method_3 = curvature_method_3 * pixel_spacing_cm_X
  ) %>%
  # Remove records with missing calibration data
  na.omit()

# Apply same calibration to DICOM-only measurements
dicom_measurements_combined <- dicom_measurements_combined %>%
  mutate(
    curvature_method_2 = curvature_method_2 * pixel_spacing_cm_X,
    cumulative1 = cumulative1 * pixel_spacing_cm_X,
    cumulative2 = cumulative2 * pixel_spacing_cm_X,
    curvature_method_3 = curvature_method_3 * pixel_spacing_cm_X
  )

# ==============================================================================
# DUPLICATE ANALYSIS - ALL MEASUREMENTS
# ==============================================================================

# Identify patients with multiple measurements across all datasets
duplicate_eids_all <- names(table(all_measurements$eid))[
  table(all_measurements$eid) > 1
]

# Create comprehensive duplicate comparison
duplicate_comparison_all <- all_measurements[
  all_measurements$eid %in% duplicate_eids_all, 
] %>%
  group_by(eid) %>%
  summarise(
    measurement_combined = list(cumulative1),
    image_name_combined = list(image_name)
  ) %>%
  # Parse combined measurements
  separate(measurement_combined, sep = ",", into = c("first", "second")) %>%
  separate(image_name_combined, sep = ",", into = c("first_image", "second_image")) %>%
  mutate(
    first = as.numeric(str_remove_all(first, "c\\(")),
    second = as.numeric(str_remove_all(second, "\\)")),
    difference = abs(first - second)
  ) %>%
  na.omit()

# Assess correlation and create visualization
correlation_all_duplicates <- cor(duplicate_comparison_all$first, 
                                  duplicate_comparison_all$second, 
                                  use = "complete.obs")
print(paste("Correlation between all duplicate measurements:", 
            round(correlation_all_duplicates, 3)))

duplicate_plot_all <- ggplot(duplicate_comparison_all, 
                             aes(x = first, y = second, color = difference)) +
  geom_point() +
  labs(
    title = "Duplicate Measurement Agreement - All Datasets",
    x = "First Measurement (cumulative1)",
    y = "Second Measurement (cumulative1)",
    color = "Absolute Difference"
  ) +
  theme_minimal()

print(duplicate_plot_all)

# ==============================================================================
# FINAL DATA CLEANING
# ==============================================================================

# Define additional problematic images for removal
additional_problematic_images <- c(
  "1.2.840.113619.2.110.212401.20230510153222.1.8.12.1",
  "1.2.840.113619.2.110.212038.20220811144246.1.12.12.1"
)

# Remove problematic images from final dataset
all_measurements <- all_measurements %>%
  filter(!str_detect(image_name, paste(additional_problematic_images, collapse = "|")))

# ==============================================================================
# LONGITUDINAL DATA ANALYSIS
# ==============================================================================

# Identify individuals with multiple measurements for longitudinal analysis
duplicate_individuals <- all_measurements %>%
  group_by(eid) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  # Add measurement visit identifier
  group_by(eid) %>%
  mutate(measurement_id = row_number()) %>%
  ungroup() %>%
  na.omit()

# Create wide format for longitudinal analysis
longitudinal_measurements_wide <- duplicate_individuals %>%
  pivot_wider(
    names_from = measurement_id,
    values_from = c(curvature_method_2, cumulative1, curvature_method_3, 
                    cumulative2, pixel_spacing_cm_X),
    names_sep = "_visit_"
  )

# Alternative wide format (all variables)
all_measurements_wide <- duplicate_individuals %>%
  pivot_wider(
    names_from = measurement_id,
    values_from = c(curvature_method_2, cumulative1, curvature_method_3, 
                    cumulative2, pixel_spacing_cm_X)
  )

# ==============================================================================
# GWAS DATASET PREPARATION AND ENHANCEMENT
# ==============================================================================

# Enhance GWAS dataset with DICOM metadata
gwas_dicom_enhanced <- gwas_dicom_8_19 %>%
  left_join(dicom_img_metadata, by = "image_name")

# Extract unique EIDs for analysis
unique_gwas_eids <- unique(gwas_dicom_enhanced$eid)
gwas_eids_df <- data.frame(gwas_eids = unique_gwas_eids)

# Create unique image dataset
unique_gwas_images <- data.frame(image_name = unique(gwas_dicom_enhanced$image_path)) %>%
  left_join(gwas_dicom_enhanced, by = c("image_name" = "image_path"))

# ==============================================================================
# CLEAN GWAS DATASET CREATION
# ==============================================================================

# Remove duplicates and create clean GWAS dataset
gwas_dicom_clean <- gwas_dicom_12_9 %>%
  group_by(image_name) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  distinct(image_name, .keep_all = TRUE)

# Create comprehensive EID dataset
eid_comprehensive <- Patient_EID %>%
  full_join(dicom_img_metadata, by = "image_name") %>%
  group_by(image_name) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  distinct(image_name, .keep_all = TRUE)

# Enhance GWAS dataset with patient information
gwas_dicom_final <- gwas_dicom_8_19 %>%
  left_join(eid_comprehensive, by = "image_name") %>%
  mutate(eid = as.numeric(eid)) %>%
  # Add sex information
  left_join(fid31_sex, by = "eid") %>%
  rename(sex = `31-0.0`) %>%
  # Add body region classification
  left_join(patient_body_region_classify_all_in_, by = "image_name")

# ==============================================================================
# DUPLICATE HANDLING IN GWAS DATASET
# ==============================================================================

# Identify and process duplicates in GWAS dataset
gwas_duplicates <- gwas_dicom_final %>%
  group_by(eid) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  select(eid, 1:14) %>%
  group_by(eid) %>%
  mutate(measurement_id = row_number()) %>%
  ungroup() %>%
  filter(!is.na(eid))

# Create wide format for duplicates
gwas_duplicates_wide <- gwas_duplicates %>%
  pivot_wider(
    names_from = measurement_id,
    names_sep = "_",
    values_from = -eid
  )

# Prepare body region classification data
patient_body_region_classify_clean <- patient_body_region_classify_all_in_ %>%
  select(1:2) %>%
  rename(image_name = file_name)

# Export duplicates for external analysis
write.csv(gwas_duplicates_wide, "gwas_duplicates_analysis.csv", row.names = FALSE)

# Create unique GWAS dataset (one record per individual)
gwas_dicom_unique <- gwas_dicom_final %>%
  distinct(eid, .keep_all = TRUE) %>%
  # Add age information
  left_join(fid21003, by = "eid")

# Add scoliosis diagnosis information
scoliosis_diagnoses <- binary_ICD_011924 %>%
  rename(eid = Patient.EID) %>%
  select(eid, M41)

gwas_dicom_with_diagnoses <- gwas_dicom_unique %>%
  left_join(scoliosis_diagnoses, by = "eid")

# ==============================================================================
# GWAS RESULTS PROCESSING
# ==============================================================================

# Process GWAS results for cumulative1 rank transform
gwas_cumulative1_results <- X70K_GWAS_090424_cumulative1ranktransform_glm_linear %>%
  # Standardize column names
  rename(
    chr = CHROM,
    pos = POS,
    rsid = ID,
    A2 = REF,
    beta = BETA,
    N = OBS_CT
  ) %>%
  # Remove unnecessary columns
  select(-ALT, -TEST, -L95, -U95, -T_STAT, -P) %>%
  # Reorder columns
  select(chr, pos, rsid, A1, A2, everything())

# Create SNP identifier for POPCORN analysis
ranktransform_popcorn <- ranktransform_popcorn %>%
  mutate(SNP = paste(chr, pos, sep = ":")) %>%
  select(chr, pos, SNP, everything()) %>%
  select(-rsid)

# Standardize column names for POPCORN analysis
ranktransform_popcorn2 <- ranktransform_popcorn2 %>%
  rename(a1 = A1, a2 = A2) %>%
  select(-chr, -pos)

# ==============================================================================
# SCOLIOSIS GWAS RESULTS PROCESSING
# ==============================================================================

# Process scoliosis GWAS results
scoliosis_gwas_results <- scoliosis_step2_M41_no_imaging_pvalue %>%
  rename(
    chr = CHROM,
    pos = GENPOS,
    a1 = ALLELE0,
    a2 = ALLELE1,
    beta = BETA
  ) %>%
  select(-ID, -A1FREQ, -TEST, -CHISQ, -LOG10P, -PVALUE) %>%
  mutate(SNP = paste(chr, pos, sep = ":")) %>%
  select(SNP, everything()) %>%
  select(-chr, -pos)

# ==============================================================================
# JAPANESE GWAS DATA PROCESSING
# ==============================================================================

# Process Japanese GWAS meta-analysis results
japanese_gwas_results <- AIS_ImputationMETA.RSQR03_MAF005 %>%
  rename(
    chr = CHR,
    pos = BP,
    `p-value` = P
  ) %>%
  select(-P.R., -OR.R., -Q, -I) %>%
  mutate(N = 79211) %>%
  # Move p-value to last column
  select(-`p-value`, everything(), `p-value`)

# Process Japanese POPCORN data
japanese_popcorn <- japanese_popcorn %>%
  rename(a1 = A1, a2 = A2, `p-value` = p.value) %>%
  select(-chr, -pos)

# Calculate beta and standard error from OR and p-value for Japanese data
japanese_popcorn2 <- japanese_popcorn2 %>%
  mutate(
    beta = log(OR),
    z_value = qnorm(1 - p.value / 2),
    SE = abs(beta / z_value)
  ) %>%
  select(SNP, a1, a2, N, beta, SE)

# ==============================================================================
# EUROPEAN-EAST ASIAN ANALYSIS
# ==============================================================================

# Process European-East Asian combined genetic data
eur_eas_combined <- EUR_EAS_all_gen_imp %>%
  mutate(SNP = paste(V1, V2, sep = ":")) %>%
  select(V1, V2, SNP, everything()) %>%
  select(-V3)

# ==============================================================================
# ADDITIONAL PHENOTYPE INTEGRATION
# ==============================================================================

# Add bone mineral density and body composition measurements
gwas_dicom_enhanced_phenotypes <- gwas_dicom_12_93 %>%
  left_join(fid22433, by = "eid") %>%    # Abdominal weight-muscle ratio
  left_join(fid23234, by = "eid") %>%    # Spine BMD
  left_join(fid23285, by = "eid") %>%    # Trunk lean mass
  left_join(fid47, by = "eid") %>%       # Right grip strength
  # Rename variables for clarity
  rename(
    spine_bmd = `23234-2.0`,
    trunk_lean_mass = `23285-2.0`,
    abdominal_weight_muscle_ratio = `22433-2.0`,
    right_grip_strength = `47-0.0`,
    sex = sex
  ) %>%
  # Convert sex to factor
  mutate(sex = as.factor(sex))

# ==============================================================================
# SUMMARY STATISTICS AND VALIDATION
# ==============================================================================

# Print comprehensive summary information
cat("COMPREHENSIVE DATASET SUMMARY\n")
cat("=============================\n")
cat("Original interpretations:", nrow(frontera_method2_interpretations_backup), "\n")
cat("New DICOM measurements:", nrow(dicom_measurements_combined), "\n")
cat("Combined measurements:", nrow(all_measurements), "\n")
cat("Unique individuals:", length(unique(all_measurements$eid)), "\n")
cat("Individuals with multiple measurements:", length(duplicate_eids_all), "\n")
cat("GWAS dataset size:", nrow(gwas_dicom_final), "\n")
cat("Unique GWAS individuals:", nrow(gwas_dicom_unique), "\n")
cat("Longitudinal measurements available:", nrow(duplicate_individuals), "\n")

# Correlation summary
if(nrow(duplicate_comparison_all) > 0) {
  cat("Overall duplicate measurement correlation:", 
      round(correlation_all_duplicates, 3), "\n")
}

cat("\nData processing completed successfully!\n")