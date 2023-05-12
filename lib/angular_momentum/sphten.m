
function V = sphten(v)

dim = find(size(v)==3);
switch dim
    case 1
        V = cat(1,-(v(1,:,:,:)+1i*v(2,:,:,:))/sqrt(2),v(3,:,:,:),(v(1,:,:,:)-1i*v(2,:,:,:))/sqrt(2));
        V = flip(V,1);
    case 2
        V = cat(2,-(v(:,1,:,:)+1i*v(:,2,:,:))/sqrt(2),v(:,3,:,:),(v(:,1,:,:)-1i*v(:,2,:,:))/sqrt(2));
        V = flip(V,2);
    case 3
        V = cat(3,-(v(:,:,1,:)+1i*v(:,:,2,:))/sqrt(2),v(:,:,3,:),(v(:,:,1,:)-1i*v(:,:,2,:))/sqrt(2));
        V = flip(V,3);
    case 4
        V = cat(4,-(v(:,:,:,1)+1i*v(:,:,:,2))/sqrt(2),v(:,:,:,3),(v(:,:,:,1)-1i*v(:,:,:,2))/sqrt(2));
        V = flip(V,4);
end

end