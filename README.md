# This is the source code to generate the results for the research: "An optimal power allocation for multi-cell downlink NOMA system" [1].

# Please follow the instructions below in order to regenerate the results:

First of all, you should generate the channels by executing the file "generate_channels.m", which saves the channels to the folders: "channels_for_NU" and "channels_for_powers", which are the channels needed to plot the sum-rate VS NU/NC and the sum-rate VS P(max), respectively, as well as the convergence plot.

After generating the channels, execute the files: "execute_WMMSE_plot_conv.m", "execute_WMMSE_plot_RvsNU.m", and "execute_WMMSE_plot_RvsPm.m" in order to generate the convergence plot, the sum-rate VS NU/NC plot, and the sum-rate VS P(max) plot, respectively. Each of these three files execute the "execute_WMMSE.m" file, which in turn executes the WMMSE algorithm and saves the data to be used later inside the folders "WMMSE_for_conv", "WMMSE_for_NU", and "WMMSE_for_powers", respectively.

To not consume time, you can download all the pre-saved data described above from the following link:
https://bit.ly/3BgvXeL



And Thank You for your interest in our research.

[1]    Raed S. H. AL-Musawi, Mujtaba A. Kareem, Ali Z. Al-Saygh, and Yahya J. Harbi, "An optimal power allocation for multi-cell downlink NOMA system," _Draft._
