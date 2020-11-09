function [ImSegmented,ImMask] = segmentObjects(Im,ImBgRemoved,ImGaussFiltered,thresh)
%parameter thresh to vary

ImMask = zeros(size(ImBgRemoved),class(ImBgRemoved));
ImMask(ImBgRemoved>(ImGaussFiltered/thresh))=1;
ImMask(ImBgRemoved<=(ImGaussFiltered/thresh))=0;

ImSegmented = ImMask.*Im;

end