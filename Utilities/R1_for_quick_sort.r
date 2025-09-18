print("Runoob")

quick_sort <- function(arr) {
  if (length(arr) <= 1) {
    return(arr)
  } else {
    pivot <- arr[length(arr) %/% 2 + 1]
    left <- arr[arr < pivot]
    middle <- arr[arr == pivot]
    right <- arr[arr > pivot]
    return(c(quick_sort(left), middle, quick_sort(right)))
  }
}

# Example usage:
arr <- c(3, 6, 8, 10, 1, 2, 1)
sorted_arr <- quick_sort(arr)
print(paste("Sorted array:", sorted_arr))
