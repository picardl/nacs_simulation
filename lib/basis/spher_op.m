function out = spher_op_v2(r,c,p,J,M)



out = (-1).^(r.(J)-r.(M)).*w3j(r.(J),-r.(M),1,p,c.(J),c.(M)).*sqrt(r.(J).*(r.(J)+1).*(2*r.(J)+1)).*(r.(J)==c.(J));

end