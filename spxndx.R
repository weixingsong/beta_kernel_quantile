# Install necessary packages (if not already installed)

# Load necessary libraries

# Load packages

library(patchwork)
library(data.table)
library(ggplot2)
library(splines)
library(transport)
# Load the package
library(quantmod)
# Get historical data for S&P 500 (^GSPC) from Yahoo Finance
getSymbols("^GSPC", src = "yahoo", from = "2023-01-01", to = "2024-12-31")
# View the first few rows
#head(GSPC)
# Convert to data frame
spx_data <- data.frame(Date = index(GSPC), coredata(GSPC))

# Download historical data from Yahoo Finance
getSymbols("^NDX", src = "yahoo", from = "2023-01-01", to = "2024-12-31")

# Convert to data frame
ndx_data <- data.frame(Date = index(NDX), coredata(NDX))


# Function to compute daily returns
calculate_daily_returns <- function(prices) {
  return(diff(log(prices), lag=1))
}

# Load SPX and NDX historical price data from CSV files
# Ensure the CSV files contain "Date" and "Close" columns
#spx_data <- fread("SPX_data.csv")
#ndx_data <- fread("NDX_data.csv")

# Convert Date column to Date type
spx_data$Date <- as.Date(spx_data$Date)
ndx_data$Date <- as.Date(ndx_data$Date)

# Calculate daily log returns
spx_returns <- calculate_daily_returns(spx_data$GSPC.Close)
ndx_returns <- calculate_daily_returns(ndx_data$NDX.Close)

# Compute empirical quantile functions
quantiles <- seq(0.0001, 0.9999, length.out = length(spx_returns))
spx_empirical <- quantile(spx_returns, probs = quantiles)
ndx_empirical <- quantile(ndx_returns, probs = quantiles)

# Fit smoothed quantile functions using cubic spline interpolation
#spx_spline <- smooth.spline(quantiles, spx_empirical_qf, spar = 0.7)
#ndx_spline <- smooth.spline(quantiles, ndx_empirical_qf, spar = 0.7)

# fit smoothed quantile functions using beta, normal kernel

cv=1

xord = sort(spx_returns)
yord = sort(ndx_returns)
n=length(spx_returns)
iseq = seq(n)/n
jseq = (seq(n)-1)/n         

Qbtspx=Qntspx=Qhtspx=QEptspx=rep(0,length(quantiles))
Qbtndx=Qntndx=Qhtndx=QEptndx=rep(0,length(quantiles))

hb=hv=cv*n^(-2/3)
hn=he=cv*n^(-1/3)
jj=1
for(pv in quantiles)
{
  Bv=pbeta(iseq,pv/hb+1,(1-pv)/hb+1)-
  pbeta(jseq,pv/hb+1,(1-pv)/hb+1)
  Nv=pnorm(iseq,pv,hn) - pnorm(jseq,pv,hn)
  Hv = pbeta(iseq,pv/hv,(1-pv)/hv)-
     pbeta(jseq,pv/hv,(1-pv)/hv)
  Sv=rep(0,n)
  for(kk in seq(n))
  {
    a=pv-hn
    b=pv+hn
    c=(kk-1)/n
    d=kk/n
    
    lend=max(a,c)
    rend=min(b,d)
    
    
    if(lend>=rend)
    {Sv[kk]=0}
    if(lend<rend)
    {
      funu=function(u)
      {
        0.75*(1-((u-pv)/hn)^2)/hn
      }
      Sv[kk]=integrate(funu,lend,rend)$value
    }
  }
 
  Qbtspx[jj] = sum(xord*Bv)
  Qntspx[jj] = sum(xord*Nv)
  Qhtspx[jj] = sum(xord*Hv)
  QEptspx[jj]= sum(xord*Sv)
  
  Qbtndx[jj] = sum(yord*Bv)
  Qntndx[jj] = sum(yord*Nv)
  Qhtndx[jj] = sum(yord*Hv)
  QEptndx[jj]= sum(yord*Sv)
  
  jj=jj+1
}  



# Compute Wasserstein-1 distance using numerical integration
W1_empirical <- mean(abs(sort(spx_empirical)-sort(ndx_empirical)))
W1_smoothed_b1 <- mean(abs(sort(Qbtspx)-sort(Qbtndx)))
W1_smoothed_b2 <- mean(abs(sort(Qhtspx)-sort(Qhtndx)))
W1_smoothed_ep <- mean(abs(sort(QEptspx)-sort(QEptndx)))
W1_smoothed_n <- mean(abs(sort(Qntspx)-sort(Qntndx)))


# Plot empirical and smoothed quantile functions
df_plot <- data.frame(
  Quantile = quantiles,
  SPX_Empirical = spx_empirical,
  NDX_Empirical = ndx_empirical,
  SPX_Smoothedb1 = Qbtspx,
  NDX_Smoothedb1 = Qbtndx,
  SPX_Smoothedb2 = Qhtspx,
  NDX_Smoothedb2 = Qhtndx,
  SPX_Smoothedep = QEptspx,
  NDX_Smoothedep = QEptndx,
  SPX_Smoothedn = Qntspx,
  NDX_Smoothedn = Qntndx
)

# First Plot: Smoothed B1
p1 <- ggplot(df_plot, aes(x = Quantile)) +
  geom_line(aes(y = SPX_Empirical), linetype = "dotted", size = 1, color = "darkgrey") +
  geom_line(aes(y = NDX_Empirical), linetype = "dotted", size = 1, color = "black") +
  geom_line(aes(y = SPX_Smoothedb1), size = 0.5, color = "darkgrey") +
  geom_line(aes(y = NDX_Smoothedb1), size = 0.5, color = "black") +
  labs(x = "", y = "Return") +
  theme_classic() +
  theme(legend.position = "none")

# Second Plot: Smoothed B2
p2 <- ggplot(df_plot, aes(x = Quantile)) +
  geom_line(aes(y = SPX_Empirical), linetype = "dotted", size = 1, color = "darkgrey") +
  geom_line(aes(y = NDX_Empirical), linetype = "dotted", size = 1, color = "black") +
  geom_line(aes(y = SPX_Smoothedb2), size = 0.5, color = "darkgrey") +
  geom_line(aes(y = NDX_Smoothedb2), size = 0.5, color = "black") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(legend.position = "none")

# Third Plot: Smoothed EP
p3 <- ggplot(df_plot, aes(x = Quantile)) +
  geom_line(aes(y = SPX_Empirical), linetype = "dotted", size = 1, color = "darkgrey") +
  geom_line(aes(y = NDX_Empirical), linetype = "dotted", size = 1, color = "black") +
  geom_line(aes(y = SPX_Smoothedep), size = 0.5, color = "darkgrey") +
  geom_line(aes(y = NDX_Smoothedep), size = 0.5, color = "black") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(legend.position = "none")

# Third Plot: Smoothed N
p4 <- ggplot(df_plot, aes(x = Quantile)) +
  geom_line(aes(y = SPX_Empirical), linetype = "dotted", size = 1, color = "darkgrey") +
  geom_line(aes(y = NDX_Empirical), linetype = "dotted", size = 1, color = "black") +
  geom_line(aes(y = SPX_Smoothedn), size = 0.5, color = "darkgrey") +
  geom_line(aes(y = NDX_Smoothedn), size = 0.5, color = "black") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(legend.position = "none")

# Combine the three plots into one figure
(p1 | p2 | p3 | p4)  # Arranges the plots horizontally


c(W1_empirical,W1_smoothed_b1,W1_smoothed_b2,W1_smoothed_ep, W1_smoothed_n)

