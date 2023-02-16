
% Load some filenames etc.
c = metadata;

% Number of timesteps per day that we analyse (up to 144)
nt = 144;

% Number of daily files that we analyse (up to 14)
nf = 5;

% Box over which we perform the analysis
x = [145 155];
y1 = [-40 -30];

% Get the latitude and longitude and land mask
orog_name = [c.output_dir '/' c.orog];

lat = ncread(orog_name,'latitude');
nlat = length(lat);


lon = ncread(orog_name,'longitude');
nlon = length(lon);

orog = ncread(orog_name,'orog',[1 1 1],[nlon nlat 1]);
land = orog>0;

% Create 2-D lat and lon matrices
[yy,xx] = meshgrid(lat,lon);

% I1 is an indicator that tells us whether we are in the analysis region
I1 = xx>x(1) & xx<=x(2) & yy>y1(1) & yy<=y1(2);

% Now lets calculate the precip distribution
percentiles = [0:0.01:.99 0.992 0.994 0.996 0.998 0.999 0.9992 0.9994 0.9996 0.9998 0.9999];
prec_dist_land = zeros(1,length(percentiles));
prec_dist_ocean = zeros(1,length(percentiles));

prec_mean_land = 0;
prec_mean_ocean = 0;
nsnap = 0;

for k = 1:nf
disp(['file ' num2str(k)])


   file_name = [c.output_dir '/' c.precip_files{k}];
   for i = 1:nt
      disp(num2str(i))
      prec = ncread(file_name,'accum_ls_prcp',[1 1 i],[nlon nlat 1]);


      % Calculate the percentiles for this snapshot
      % this script calculates averages of the percentiles for each snapshot
      % It would be better to calculate the percentiles over all snapshots
      prec_dist_ocean = prec_dist_ocean + quantile(prec(land<0.5 & I1),percentiles);
      prec_dist_land = prec_dist_land + quantile(prec(land>0.5 & I1),percentiles);

      prec_mean_land = prec_mean_land + mean(prec(land<0.5 & I1));
      prec_mean_ocean = prec_mean_ocean + mean(prec(land>0.5 & I1));
 
      nsnap = nsnap + 1;
  
   end

end

prec_dist_ocean = prec_dist_ocean./nsnap;
prec_dist_land = prec_dist_land./nsnap;
   

% Now plot the precipitation distribution
fig.bfig(15,10)

plot(prec_dist_ocean,1-percentiles,'o-')
hold on
plot(prec_dist_land,1-percentiles,'o-')

set(gca,'ytick',[0.0001 0.001 0.01 0.1],'yticklabel',{'99.99th' '99.9th' '99th' '90th'})

xlabel('Precipitation rate')
ylabel('percentile') 

set(gca,'yscale','log')

save('precip_dist.mat','precip_dist_land','precip_dist_ocean','percentiles','lat','lon','x','y')


print -dpdf ./precip_dist.pdf

% Also make a basic plot of the analysis region
fig.bfig(15,10)
[xc,yc] = earth.coasts;
pcolor(lon,lat,I1'.*1);
shading flat;
hold on
plot(xc,yc);
xlabel('longitude (deg)')
ylabel('latitudfe (deg)')


print -dpdf ./map.pdf



