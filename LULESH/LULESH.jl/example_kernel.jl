import JACC

function axpy(i, alpha, x, y)
    if i <= length(x)
        @inbounds x[i] += alpha * y[i]
    end
end

N = 10
# Generate random vectors x and y of length N for the interval [0, 100]
x = round.(rand(Float32, N) * 100)
y = round.(rand(Float32, N) * 100)
alpha = 2.5

x_d = JACC.Array(x)
y_d = JACC.Array(y)
# Execute 1000 times
for i in 1:1000000
    JACC.parallel_for(N, axpy, alpha, x_d, y_d)
end
# To verify the result, you can print the updated x_d array
println("Updated x_d: ", Array(x_d))