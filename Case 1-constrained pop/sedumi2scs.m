function [data,K] = sedumi2scs(At,b,C,K)

%% SEDUMI2SCS.m Convert sedumi input to SCS

% Get non SDP variables
nonSDP = 0;
if isfield(K,'f') && ~isempty(K.f); nonSDP = nonSDP + K.f; else K.f = 0; end
if isfield(K,'l') && ~isempty(K.l); nonSDP = nonSDP + K.l; else K.l = 0;end
if isfield(K,'q') && ~isempty(K.q); nonSDP = nonSDP + sum(K.q); else K.q = 0;end
if isfield(K,'r') && ~isempty(K.r); nonSDP = nonSDP + sum(K.r); end

% Turn At for SDP variables from vec to svec
if isfield(K,'s') && ~isempty(K.s)
    At0 = At(1:nonSDP,:);
    Ats = At(nonSDP+1:end,:);
    C0 = C(1:nonSDP);
    Cs = C(nonSDP+1:end);
    count = 0;
    for i = 1:length(K.s)
       Q = svecTransMat(K.s(i));
       At0 = [At0; Q*Ats(count+1:count+K.s(i)^2,:)];
       C0 = [C0; Q*Cs(count+1:count+K.s(i)^2)];
       count = count+K.s(i)^2;
    end
end

% Set outputs
data.A = At0;
data.b = full(C0);
data.c = -full(b);
