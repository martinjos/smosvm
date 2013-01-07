function feats = findbestfeatures(w)

feats = 1:length(w);

for i=1:size(w, 1)
    feats = intersect(feats, find(w(i,:) > mean(w(i,:))));
end

disp(length(feats));