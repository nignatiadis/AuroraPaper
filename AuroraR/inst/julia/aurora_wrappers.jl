function auroral(Zs)
  res = StatsBase.fit(Auroral(), ReplicatedSample.(Zs))
  Dict(:mus => res.μs,
       :betas => res.βs,
       :mus_mat => res.μs_mat,
       :betas_list => res.βs_list)
end

function ebcc(Zs)
  res = StatsBase.fit(CoeyCunningham(), ReplicatedSample.(Zs))
  Dict(:mus => res.μs,
       :betas => res.βs,
       :mus_mat => res.μs_mat,
       :betas_list => res.βs_list)
end

function aurora_knn(Zs, kKNN, loocv)
  res = StatsBase.fit(AuroraKNN(kKNN=kKNN, loocv=loocv), ReplicatedSample.(Zs))
  Dict(:mus => res.μs,
       :mus_mat => res.μs_mat,
       :ks => res.kskNN)
end
