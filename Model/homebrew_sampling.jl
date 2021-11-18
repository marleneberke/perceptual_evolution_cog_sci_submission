#Sample n times from the array_to_sample_from with replacement
function homebrew_sample(array_to_sample_from::Array, n::Int)
    array_to_return = []
    for i=1:n
        push!(array_to_return, rand(array_to_sample_from))
    end
    return array_to_return
end

#Sample n times from the array_to_sample_from without replacement
function sample_wo_repl!(A::Array, n::Int)
    sample = []
    for i in 1:n
        push!(sample, splice!(A, rand(eachindex(A))))
    end
    return sample
end

#Sample n times from the array_to_sample_from with replacement. Weight each element
#in array_to_sample_from with a weight
function homebrew_sample(array_to_sample_from::Array{Int64}, weights::Array, n::Int)
    array_to_return = Array{Int64}(undef, n)

    weights[weights.<0] .= 0 #to deal with negative fitnesses when there are penalties

    if sum(weights)==0 #if all the payoffs were 0, make all the fitnesses equal
        weights = repeat([1], n_players)
    end

    weights = weights/sum(weights) #normalize just in case

    # for i=1:n
    #     r = rand() #produce a random number between 0 and 1
    #     w = 1
    #     while r > sum(weights[1:w])
    #         w = w+1
    #     end
    #     push!(array_to_return, array_to_sample_from[w])
    # end
    # return array_to_return
    #problem was it was too slow

    dist = Categorical(weights)
    for i=1:n
        array_to_return[i] = array_to_sample_from[rand(dist)]
    end
    return array_to_return
end
