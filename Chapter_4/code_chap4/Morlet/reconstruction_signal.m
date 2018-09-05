function signal_cont=reconstruction_signal(signal_discr,t_step,t_discr_first,t_cont)

    % Function used to reconstruct, at any time, the signal from its samples
    % with the "sinc" function (see the Shannon theory on continuous signals with
    % finite frequency bandwidth). 
    
    % signal_discr: the signal sampled at t_step
    % t_step: time step of the sampled signal
    % t_discr_first: first time of the sampled signal
    % t_cont: the times at which we want the reconstructed signal
    
    mysize=size(signal_discr);
    if mysize(1)==1
        signal_discr=transpose(signal_discr);
    end
    N=numel(signal_discr);
    myfun=@(t) sinc((t-t_discr_first)./t_step-(0:N-1))*signal_discr;
    N_cont=numel(t_cont);
    signal_cont=zeros(N_cont,1);
    for k=1:N_cont
        signal_cont(k)=myfun(t_cont(k));
    end

end