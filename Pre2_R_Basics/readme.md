This course references content from *Harvard Stat 115/215 Lab1 R Basics*, as a prerequisite for quickly getting started with Bioinformatics.

The corresponding Jupyter Notebook file is also available for direct access.

## Table of Contents
  - [1 Install Relevant R Packages](#1-install-relevant-r-packages)
  - [2 Data Types and Data Structures in R](#2-data-types-and-data-structures-in-r)
    - [Basic Data Types in R:](#basic-data-types-in-r)
    - [Common Data Structures in R:](#common-data-structures-in-r)
  - [3 Vector in R](#3-vector-in-r)
  - [4 Matrix in R](#4-matrix-in-r)
  - [5 List in R](#5-list-in-r)
  - [6 Data Frames in R](#6-data-frames-in-r)
  - [7 Data Manipulation with `dplyr`](#7-data-manipulation-with-dplyr)
    - [7.1 filter (Row Selection by Condition)](#71-filter-row-selection-by-condition)
    - [7.2 select (Column Selection)](#72-select-column-selection)
    - [7.3 mutate (Adding New Columns)](#73-mutate-adding-new-columns)
    - [7.4 arrange (Sorting Data)](#74-arrange-sorting-data)
    - [7.5 summarize (Aggregation)](#75-summarize-aggregation)
    - [7.6 group_by (Grouping Data)](#76-group_by-grouping-data)
  - [8 Workflow for Data Processing](#8-workflow-for-data-processing)
  - [9 The Pipe Operator `%>%`](#9-the-pipe-operator-)
    - [Using Pipe Operator with `group_by`](#using-pipe-operator-with-group_by)
    - [Solution to the Exercise](#solution-to-the-exercise)
  - [10 Visualization with `ggplot2`](#10-visualization-with-ggplot2)
    - [Example 1](#example-1)
    - [Code Explanation](#code-explanation)
    - [Example 2](#example-2)

## 1 Install Relevant R Packages

In a previous article (Pre1), we created an R environment using Anaconda, installed the necessary packages, and configured Jupyter Notebook to support R.

This section requires three R packages: `ggplot2`, `dplyr`, and `nycflights13`.

Since `ggplot2` and `dplyr` are already included in `r-essentials`, we only need to install `nycflights13`.

```bash
conda activate <environment_name>
conda install -c conda-forge r-nycflights13
```

**Note:** The `nycflights13` package provides a well-known dataset collection designed for learning and teaching data analysis. It records all **flights** from the three major New York airports (JFK, LGA, EWR) in 2013. It is often used with `dplyr` and `ggplot2` for teaching data cleaning and exploratory data analysis.

## 2 Data Types and Data Structures in R

R’s basic data types and structures.

### Basic Data Types in R:

- **Character** (`"A"`)
- **Numeric** (`1`)
- **Integer** (`1L`)
- **Logical** (`TRUE`, `T`)
- **Complex** (`1+1i`)

### Common Data Structures in R:

- **Atomic Vectors**
- **Lists**
- **Matrices** (Matrix)
- **Data Frames**
- **Factors**

> **Note:**
>
> ① **What is *Atomic*:**
>
> - A vector or collection is considered atomic if all elements have the same data type.
>
> ② **What is *Factor*:**
>
> - **Factors** in R are a special data structure primarily used to represent **categorical data**.
> - Factors store distinct categories as **discrete levels**, using integer values to represent them. This optimizes memory usage and enhances processing efficiency.

## 3 Vector in R

R natively supports vectorized operations, which differs from **List** in *Python*.

```R
# Creating a vector
# In R, the assignment operator is `<-`
x <- c(1, 2, 3, 4, 5)
print(x[2]) # Accessing elements by index (indexing starts from 1)

# Performing operations on vector elements
x^2
sqrt(x)

# Logical vectors
print(x[x < 3]) # Using logical indexing

# Initializing an empty numeric vector
vector("numeric", 5)
```

![Image](https://github.com/user-attachments/assets/1d840d74-7608-4aba-99d5-b53d88bcad11)

## 4 Matrix in R

A matrix is an atomic vector with dimension attributes.

```R
# Creating a 2x4 matrix
y <- matrix(1:8, nrow = 2, ncol = 4, byrow = FALSE)
y # View the matrix
str(y) # Check matrix structure

# Accessing matrix elements
y[1, 2] # First row, second column
y[, 2]  # All elements in the second column
dim(y)  # Matrix dimensions
y %*% t(y) # Matrix multiplication
```

![Image](https://github.com/user-attachments/assets/e962f903-d91e-4667-8382-af9f30cf5cc3)

```R
dim(y) <- NULL  # Removing dimension attributes flattens the matrix into a vector
y
```

![Image](https://github.com/user-attachments/assets/51fd12f6-bf1d-41a3-b0f8-ed29d1cff766)

## 5 List in R

Lists are generalized vectors that can contain elements of different types.

```R
# Creating a list containing a vector, matrix, and another list
list_data <- list(c("Jan", "Feb", "Mar"),
                  matrix(c(3, 9, 5, 1, -2, 8), nrow = 2),
                  list("green", 12.3))

# Assigning names to list elements
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")

# Displaying the list structure
str(list_data)

# Displaying the list
list_data
```

![Image](https://github.com/user-attachments/assets/2209ce8b-e74c-4eae-9405-9fb85e78eb90)

## 6 Data Frames in R

A **data frame** is a list of equal-length R vectors.

```R
data(mtcars) # Loading a well-known dataset (from `nycflights13`)
str(mtcars)  # Checking the structure (a list of numeric vectors)
head(mtcars) # Displaying the first six rows
```

![Image](https://github.com/user-attachments/assets/33f6a2e9-6b6c-4f5b-8d5d-315dc9b9e938)

```R
mtcars[1, 1] # Accessing the value at the first row, first column
head(mtcars[1]) # Extracting the first column as a data frame
head(mtcars[[1]]) # Extracting the first column as a vector
sapply(mtcars, sum) # Computing the sum for each column
```

![Image](https://github.com/user-attachments/assets/5929c0d0-7399-4605-bd7e-924934d88c67)

## 7 Data Manipulation with `dplyr`

- `**dplyr**` is an R package often considered a "language" for data manipulation within R, making data operations intuitive and straightforward.
- In data science (and computational biology), **80% of the work is data cleaning, while only 20% is actual analysis**.
- Official `dplyr` documentation: https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html

```R
library(dplyr)
library(nycflights13)

data(flights) # Load dataset
head(flights)
```

![Image](https://github.com/user-attachments/assets/efd08e3d-f789-42a1-bf64-647ba94376f2)

### 7.1 filter (Row Selection by Condition)

- Select specific rows based on conditions.

```R
# Set the number of rows displayed in HTML tables
options(repr.matrix.max.rows = 14) # Max rows displayed
options(repr.matrix.max.cols = 8)  # Max columns displayed (optional)

filter(flights, month == 1, day == 1)
```

![Image](https://github.com/user-attachments/assets/28337fec-183f-448b-bd91-7d92eacc27f9)

### 7.2 select (Column Selection)

- Select specific columns from a dataset.

```R
select(flights, year, month, day)
```

![Image](https://github.com/user-attachments/assets/ea8e985c-763b-41b3-a8c5-97aaa09d9bea)

### 7.3 mutate (Adding New Columns)

- Add new columns derived from existing ones.

```R
mutate(flights,
  gain = arr_delay - dep_delay,
  speed = distance / air_time * 60
)
```

![Image](https://github.com/user-attachments/assets/d4a433de-dd3f-4497-af3d-8d94791183af)

### 7.4 arrange (Sorting Data)

- Sort a dataframe by specific columns.

```R
arrange(flights, year, month, day)
```

![Image](https://github.com/user-attachments/assets/7028a64b-8a57-4faf-b97f-ff0ab7a6ee34)

**Descending Order with** `desc`

- Use `desc()` to sort in descending order.

```R
select(arrange(flights, desc(dep_delay)), year, month, day, dep_delay)
```

![Image](https://github.com/user-attachments/assets/bc9adeb3-8904-42c8-9980-ab3f5de3af8f)

### 7.5 summarize (Aggregation)

- Aggregate multiple values into a single summary value.

```R
summarise(flights,
  delay = mean(dep_delay, na.rm = TRUE) # Compute mean delay, ignoring missing values
)
```

![Image](https://github.com/user-attachments/assets/a5552814-e180-4a70-a462-2238a01a6424)

### 7.6 group_by (Grouping Data)

- The true power of `dplyr` lies in its **group-apply-summarize** workflow.
- Example: Group data by aircraft tail number (`tailnum`) and compute flight count, average distance, and average delay.

```R
by_tailnum <- group_by(flights, tailnum)
delay <- summarise(by_tailnum,
   count = n(),  # Count number of flights
   dist = mean(distance, na.rm = TRUE),
   arr_delay = mean(arr_delay, na.rm = TRUE))
delay <- filter(delay, count > 20, dist < 2000)
delay
```

![Image](https://github.com/user-attachments/assets/84cb1611-70de-4b4a-901d-c25c7d3b61c7)

## 8 Workflow for Data Processing

How to apply multiple data processing steps in one go? Consider the following example:

> ```R
> a1 <- group_by(flights, year, month, day)
> a2 <- select(a1, arr_delay, dep_delay)
> a3 <- summarise(a2,
>    arr = mean(arr_delay, na.rm = TRUE),
>    dep = mean(dep_delay, na.rm = TRUE))
> a4 <- filter(a3, arr > 30 | dep > 30)
> ```

Notice that several variables (such as `a1`) are created just to be used in the next step and are never used again. This is generally not a good programming practice as it leads to unnecessary variable usage.

Now consider this:

> ```R
> filter(
>    summarise(
>        select(
>           group_by(flights, year, month, day),
>           arr_delay, dep_delay
>        ),
>        arr = mean(arr_delay, na.rm = TRUE),
>        dep = mean(dep_delay, na.rm = TRUE)
>  ),
>   arr > 30 | dep > 30
> )
> ```

This code is quite difficult to read!

The `Pipe Operator` in R offers a more structured approach for data processing workflows.

## 9 The Pipe Operator `%>%`

Writing code in a deeply nested fashion, from the innermost function to the outermost, can make it difficult to read. To solve this, the `%>%` operator (**Pipe Operator**) can be used. It transforms `f(x, y)` into `x %>% f(y)`, making the code more readable.

> ```R
> flights %>%
>    group_by(year, month, day) %>%
>    select(arr_delay, dep_delay) %>%
>    summarise(
>        arr = mean(arr_delay, na.rm = TRUE),
>        dep = mean(dep_delay, na.rm = TRUE)
>    ) %>%
>  filter(arr > 30 | dep > 30)
> ```

This approach makes the sequence of operations much clearer and more natural.

**Note:** The `%>%` pipe operator is not invented by `dplyr` but comes from the `magrittr` package. It can also be used in other scenarios, such as:

> ```R
> letters %>% length() # Compute the length of the 'letters' variable
> ```

### Using Pipe Operator with `group_by`

We can use the pipe operator together with `group_by`.

For example, consider the following code:

```R
flights %>%
  filter(origin == 'EWR') %>%
  group_by(dest) %>%
  summarize(n = n()) %>%
  arrange(desc(n))
```

![Image](https://github.com/user-attachments/assets/0a7e6177-d201-4965-b40b-12ae219fe1bb)

**Exercise:**

We have seen several functions for data transformation:

1. `select`
2. `filter`
3. `mutate`
4. `arrange`
5. `summarise`

Now, apply the following sequential operations to the `colon` dataset (make sure to use `%>%`):

1. Remove the `study` column as it contains only the value 1.
2. Retain only rows where `etype` is 2 (death events instead of recurrence events).
3. Convert `time` into months (assuming each month has exactly 30 days).
4. Arrange the dataset in ascending order by `age`.
5. Group by treatment method (`rx`) and compute the mean survival time in months.

**About `colon` dataset from the `survival` package**

The `colon` dataset is a classic dataset included in the `survival` package of R, primarily used for analyzing the survival of colon cancer patients. It contains various clinical and trial-related information, enabling the exploration of survival times in relation to different factors. Key variables include treatment method (`rx`), gender (`1 = male`), age, and trial-specific indicators. The `etype` variable indicates the type of event (death or recurrence).

```R
library(survival)
# Load dataset
data(colon)
# View first six rows
head(colon)
```

![Image](https://github.com/user-attachments/assets/1fba78da-0791-44a5-8452-b60a5f3e355e)

### Solution to the Exercise

```R
colon %>%
    select(-study) %>%
    filter(etype == 2) %>%
    mutate(time = round(time / 30)) %>%
    arrange(age) %>%
    group_by(rx) %>%
    summarize(avg_months = mean(time))
```

![Image](https://github.com/user-attachments/assets/8b4a8aef-0958-40b5-b38f-1234689ca8a3)

## 10 Visualization with `ggplot2`

We can use `ggplot2` for visualization.

### Example 1

Filter the `flights` dataset and create a **boxplot** to show arrival delays by month, with colors representing different origins.

```R
library(ggplot2)
flights %>%
    filter(arr_delay <= 360) %>%
    ggplot(aes(x = factor(month), y = arr_delay, color = origin)) +
    geom_boxplot() +
    ggtitle("Arrival Delays by Month") +
    xlab("Month") +
    ylab("Arrival Delay (minutes)")
```

![Image](https://github.com/user-attachments/assets/48a06452-9059-4dcf-9408-770788e3f38b)

### Code Explanation

**1. Initialize the plot:**

- `ggplot(aes(x = factor(month), y = arr_delay, color = origin))`:
  - `x = factor(month)`: Convert `month` into a categorical variable for the X-axis.
  - `y = arr_delay`: Map `arr_delay` to the Y-axis.
  - `color = origin`: Differentiate data points by `origin` using colors.

**2. Add a boxplot:**

- `geom_boxplot()`
  - Creates a boxplot to show the distribution of arrival delays per month.
  - The boxplot displays median, quartiles, and potential outliers.

**3. Add labels:**

- `ggtitle("Arrival Delays by Month")`: Sets the chart title.
- `xlab("Month")`: Labels the X-axis as "Month".
- `ylab("Arrival Delay (minutes)")`: Labels the Y-axis as "Arrival Delay (minutes)".

### Example 2

```R
flights %>%
    filter(month == 1, arr_delay < 360, dep_delay < 360) %>%
    ggplot(aes(x = dep_delay, y = arr_delay)) +
    geom_point() +
    ggtitle("Relation between Departure and Arrival Delay") +
    xlab("Departure Delay (min)") +
    ylab("Arrival Delay (min)")
```

![Image](https://github.com/user-attachments/assets/f9e63084-239c-43c4-8cba-5c82f5c61290)
