function f = predictedPlaceboValue(t, p)
    % p.M : mesor (daily average of rhythm pg/ml)
    % p.A : amplitude of cosine (pg/ml)
    % p.phi : acrophase (time of peak)
    f = p.M * (1 + p.A * cos((t - p.phi) * (pi/12)));
end