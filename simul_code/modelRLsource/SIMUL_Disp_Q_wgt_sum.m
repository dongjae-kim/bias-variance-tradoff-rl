%% simple Q integration (with p, n=1)
close all

if(0)
    pq_list=[1 3 10];
    for pq=pq_list
        str_box=sprintf('Q-value as a function of Q_(fwd) and Q_(sarsa) with p=%02.1f',pq);
        f2_2=figure('Name',str_box);
        f_ax2_2 = axes('Parent',f2_2);
        range_minmax=[0 20];
        resolution=100;
        Q_range=[range_minmax(1):(range_minmax(2)-range_minmax(1))/(resolution-1):range_minmax(2)];
        [Q_x,Q_y] = meshgrid(Q_range, Q_range);
        Q_z=(Q_x.^pq+Q_y.^pq).^(1/pq);
        surf(f_ax2_2 ,Q_x,Q_y,Q_z)
        view(2); colorbar;
    end
end

%% integration in arbitration (given p, varying n)
pq=[5.6];
n_list=[0.3:0.1:0.7];
for n=n_list
    str_box=sprintf('Q-value as a function of Q_(fwd) and Q_(sarsa) with n=%0.1f, p=%02.1f',n,pq);
    f2_2=figure('Name',str_box);
    f_ax2_2 = axes('Parent',f2_2);
    range_minmax=[0 20];
    resolution=100;
    Q_range=[range_minmax(1):(range_minmax(2)-range_minmax(1))/(resolution-1):range_minmax(2)];
    [Q_x,Q_y] = meshgrid(Q_range, Q_range);
    Q_z=((n*Q_x).^pq+((1-n)*Q_y).^pq).^(1/pq);
    surf(f_ax2_2 ,Q_x,Q_y,Q_z)
    view(2); colorbar;
end
