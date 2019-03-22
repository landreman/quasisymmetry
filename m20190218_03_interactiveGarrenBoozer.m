function m20190218_03_interactiveGarrenBoozer()

% Absolute path to the directory containing the fortran executable:
% Be sure to end quasisymmetry_executable_directory with '/'
quasisymmetry_executable_directory = '/Users/mattland/quasisymmetry/';

% Directory in which the input and output files will be written.
% Be sure to end quasisymmetry_temp_files_directory with '/'
quasisymmetry_temp_files_directory = '/Users/mattland/quasisymmetry/matlab/';

sign_G = 1;
sign_psi = 1;


nfp = 3;
R0c1 = 0.045;
Z0s1 = 0.045;
eta_bar = 1;
sigma = 0;
I2 = 0;
minus_mu0_p2 = 0;
B2c=0;
B2s=0;
r0 = 0.1;

X1c = 0;
Y1s = 0;
Y1c = 0;
B20 = 0;
iota = 0;
elongation = 0;
X20 = 0;
X2s = 0;
X2c = 0;
Y20 = 0;
Y2s = 0;
Y2c = 0;

Ntheta = 40;
Nphi = 150;
theta = linspace(0,2*pi,Ntheta);
phi = linspace(0,2*pi,Nphi);
[phi2D, theta2D] = meshgrid(phi,theta);
sinphi = sin(phi2D);
cosphi = cos(phi2D);

first_update = true;
number_of_calls = 0;
surf_handle = 0;
slices_plot_handle = 0;
%{
R_bold = 0;
Z_bold = 0;
R=0;
Z=0;
B=0;
%}

f = figure('Visible','off','Units','pixels','Position',[0,0,1600,800]);

ax = axes('Units','pixels','Position',[300,32,950,600]);

label_margin=25;
margin = 45;
height=615;
label_left = 30;
label_width=200;
text_left = 100;
text_width=100;
button_width = 200;
button_height = 25;

eventName = 'PostSet';

nfp_min = 1;
nfp_max = 6;

%slider_nfp = uicontrol('Style','slider','Min',0.1,'Max',2,'Value',nfp,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_nfp_callback);
slider_nfp = uicontrol('Style','slider','Min',nfp_min,'Max',nfp_max,'SliderStep',[1,1]*(1/(nfp_max-nfp_min)),'Value',nfp,'Units','pixels','Position',[30,height,200,20]);
%slider_nfp = uicontrol('Style','slider','Min',1,'Max',6,'Value',nfp,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_nfp,'Value',eventName,@slider_nfp_callback);
label_nfp = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','nfp','horizontalalignment','left');
text_nfp = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_nfp_callback);

height = height - margin;
%slider_R0c1 = uicontrol('Style','slider','Min',0.1,'Max',2,'Value',R0c1,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_R0c1_callback);
slider_R0c1 = uicontrol('Style','slider','Min',-0.3,'Max',0.3,'Value',R0c1,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_R0c1,'Value',eventName,@slider_R0c1_callback);
label_R0c1 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','R0c1','horizontalalignment','left');
text_R0c1 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_R0c1_callback);

height = height - margin;
%slider_Z0s1 = uicontrol('Style','slider','Min',0.1,'Max',3,'Value',Z0s1,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_Z0s1_callback);
slider_Z0s1 = uicontrol('Style','slider','Min',-0.3,'Max',0.3,'Value',Z0s1,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_Z0s1,'Value',eventName,@slider_Z0s1_callback);
label_Z0s1 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','Z0s1','horizontalalignment','left');
text_Z0s1 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_Z0s1_callback);

height = height - margin;
%slider_eta_bar = uicontrol('Style','slider','Min',0.1,'Max',3,'Value',eta_bar,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_eta_bar_callback);
slider_eta_bar = uicontrol('Style','slider','Min',0,'Max',3,'Value',eta_bar,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_eta_bar,'Value',eventName,@slider_eta_bar_callback);
label_eta_bar = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','eta_bar','horizontalalignment','left');
text_eta_bar = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_eta_bar_callback);

height = height - margin;
%slider_sigma = uicontrol('Style','slider','Min',-3,'Max',3,'Value',sigma,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_sigma_callback);
slider_sigma = uicontrol('Style','slider','Min',-3,'Max',3,'Value',sigma,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_sigma,'Value',eventName,@slider_sigma_callback);
label_sigma = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','sigma','horizontalalignment','left');
text_sigma = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_sigma_callback);

height = height - margin;
%slider_sigma = uicontrol('Style','slider','Min',-3,'Max',3,'Value',sigma,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_sigma_callback);
slider_I2 = uicontrol('Style','slider','Min',-3,'Max',3,'Value',I2,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_I2,'Value',eventName,@slider_I2_callback);
label_I2 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','I2','horizontalalignment','left');
text_I2 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_I2_callback);

height = height - margin;
%slider_minus_mu0_p2 = uicontrol('Style','slider','Min',0,'Max',20,'Value',minus_mu0_p2,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_minus_mu0_p2_callback);
slider_minus_mu0_p2 = uicontrol('Style','slider','Min',0,'Max',10,'Value',minus_mu0_p2,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_minus_mu0_p2,'Value',eventName,@slider_minus_mu0_p2_callback);
label_minus_mu0_p2 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','- mu0 p2','horizontalalignment','left');
text_minus_mu0_p2 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_minus_mu0_p2_callback);

height = height - margin;
%slider_B2c = uicontrol('Style','slider','Min',-3,'Max',3,'Value',B2c,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_B2c_callback);
slider_B2c = uicontrol('Style','slider','Min',-20,'Max',20,'Value',B2c,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_B2c,'Value',eventName,@slider_B2c_callback);
label_B2c = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','B2c','horizontalalignment','left');
text_B2c = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_B2c_callback);

height = height - margin;
%slider_B2s = uicontrol('Style','slider','Min',-3,'Max',3,'Value',B2s,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_B2s_callback);
slider_B2s = uicontrol('Style','slider','Min',-3,'Max',3,'Value',B2s,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_B2s,'Value',eventName,@slider_B2s_callback);
label_B2s = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','B2s','horizontalalignment','left');
text_B2s = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_B2s_callback);

height = height - margin;
%slider_r0 = uicontrol('Style','slider','Min',0,'Max',0.3,'Value',r0,'Units','pixels','Position',[30,height,200,20],'Callback',@slider_r0_callback);
slider_r0 = uicontrol('Style','slider','Min',0,'Max',0.3,'Value',r0,'Units','pixels','Position',[30,height,200,20]);
addlistener(slider_r0,'Value',eventName,@slider_r0_callback);
label_r0 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','r0','horizontalalignment','left');
text_r0 = uicontrol('Style','edit','Units','pixels','Position',[text_left,height+label_margin,text_width,20],'Callback',@text_r0_callback);

height = height - margin;
finite_r_nonlinear_button = uicontrol('Style','checkbox','Units','pixels','Position',[label_left,height+label_margin,button_width,button_height],'String','Use nonlinear finite-r method','Callback',@finite_r_nonlinear_callback,'fontsize',12);

height = height - 110;
%order_r_squared_button = uicontrol('Style','checkbox','Units','pixels','Position',[label_left,height+label_margin,button_width,button_height],'String','Include O(r^2) terms','Callback',@order_r_squared_callback,'fontsize',12);
order_r_option_button_group = uibuttongroup('Units','pixels','Position',[label_left,height+label_margin,button_width,105],'SelectionChangedFcn',@order_r_option_callback);
spacing=20;
order_r_option_r1_button = uicontrol(order_r_option_button_group,'style','radiobutton','String','O(r^1)','position',[10,5+4*spacing,200,15]);
order_r_option_r2_button = uicontrol(order_r_option_button_group,'style','radiobutton','String','O(r^2)','position',[10,5+3*spacing,200,15]);
order_r_option_r3_flux_constraint_button = uicontrol(order_r_option_button_group,'style','radiobutton','String','O(r^3), flux constraint','position',[10,5+2*spacing,200,15]);
order_r_option_r3_simplified_button = uicontrol(order_r_option_button_group,'style','radiobutton','String','O(r^3), simplified','position',[10,5+1*spacing,200,15]);
order_r_option_r3_simplified_with_Z3_button = uicontrol(order_r_option_button_group,'style','radiobutton','String','O(r^3), simplified, with Z3','position',[10,5+0*spacing,200,15]);

height = height - margin/2;
label_iota = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,label_width,20],'String','Iota: ','horizontalalignment','left','fontsize',15);

height = height - margin/2;
label_elongation = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,150,20],'String','Elongation: ','horizontalalignment','left','fontsize',15);

height = height - margin/2;
label_helicity = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,200,20],'String','Helicity: ','horizontalalignment','left','fontsize',15);

%{
height = height - margin/2;
label_B20 = uicontrol('Style','text','Units','pixels','Position',[label_left,height+label_margin,150,20],'String','B20: ','horizontalalignment','left');
%}

%save_button = uicontrol('Style','pushbutton','Units','pixels','Position',[115,40,120,50],'String','Save VMEC input file','Callback',@save_button_callback);

set_text_nfp()
set_text_R0c1()
set_text_Z0s1()
set_text_eta_bar()
set_text_sigma()
set_text_I2()
set_text_minus_mu0_p2()
set_text_B2c()
set_text_B2s()
set_text_r0()


update()

set(f,'Visible','on')

    function finite_r_nonlinear_callback(source,callbackdata)
        update()
    end
    function order_r_option_callback(source,callbackdata)
        update()
    end

    function slider_nfp_callback(source,callbackdata)
        nfp = round(slider_nfp.Value);
        set_text_nfp()
        update()
    end
    function slider_R0c1_callback(source,callbackdata)
        R0c1 = slider_R0c1.Value;
        set_text_R0c1()
        update()
    end
    function slider_Z0s1_callback(source,callbackdata)
        Z0s1 = slider_Z0s1.Value;
        set_text_Z0s1()
        update()
    end
    function slider_eta_bar_callback(source,callbackdata)
        eta_bar = slider_eta_bar.Value;
        set_text_eta_bar()
        update()
    end
    function slider_sigma_callback(source,callbackdata)
        sigma = slider_sigma.Value;
        set_text_sigma()
        update()
    end
    function slider_I2_callback(source,callbackdata)
        I2 = slider_I2.Value;
        set_text_I2()
        update()
    end
    function slider_minus_mu0_p2_callback(source,callbackdata)
        minus_mu0_p2 = slider_minus_mu0_p2.Value;
        set_text_minus_mu0_p2()
        update()
    end
    function slider_B2c_callback(source,callbackdata)
        B2c = slider_B2c.Value;
        set_text_B2c()
        update()
    end
    function slider_B2s_callback(source,callbackdata)
        B2s = slider_B2s.Value;
        set_text_B2s()
        update()
    end
    function slider_r0_callback(source,callbackdata)
        r0 = slider_r0.Value;
        set_text_r0()
        update()
    end

    function text_nfp_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = round(max([min([val,slider_nfp.Max]),slider_nfp.Min]));
            nfp = val;
            slider_nfp.Value = nfp;
            % The above line does not cause the slider callback to be called.
            set_text_nfp()
            update()
        else
            set_text_nfp()
        end
    end

    function text_R0c1_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_R0c1.Max]),slider_R0c1.Min]);
            R0c1 = val;
            slider_R0c1.Value = R0c1;
            % The above line does not cause the slider callback to be called.
            set_text_R0c1()
            update()
        else
            set_text_R0c1()
        end
    end
    function text_Z0s1_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_Z0s1.Max]),slider_Z0s1.Min]);
            Z0s1 = val;
            slider_Z0s1.Value = Z0s1;
            % The above line does not cause the slider callback to be called.
            set_text_Z0s1()
            update()
        else
            set_text_Z0s1()
        end
    end
    function text_eta_bar_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_eta_bar.Max]),slider_eta_bar.Min]);
            eta_bar = val;
            slider_eta_bar.Value = eta_bar;
            % The above line does not cause the slider callback to be called.
            set_text_eta_bar()
            update()
        else
            set_text_eta_bar()
        end
    end
    function text_sigma_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_sigma.Max]),slider_sigma.Min]);
            sigma = val;
            slider_sigma.Value = sigma;
            % The above line does not cause the slider callback to be called.
            set_text_sigma()
            update()
        else
            set_text_sigma()
        end
    end
    function text_I2_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_I2.Max]),slider_I2.Min]);
            I2 = val;
            slider_I2.Value = I2;
            % The above line does not cause the slider callback to be called.
            set_text_I2()
            update()
        else
            set_text_I2()
        end
    end
    function text_minus_mu0_p2_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_minus_mu0_p2.Max]),slider_minus_mu0_p2.Min]);
            minus_mu0_p2 = val;
            slider_minus_mu0_p2.Value = minus_mu0_p2;
            % The above line does not cause the slider callback to be called.
            set_text_minus_mu0_p2()
            update()
        else
            set_text_minus_mu0_p2()
        end
    end
    function text_B2c_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_B2c.Max]),slider_B2c.Min]);
            B2c = val;
            slider_B2c.Value = B2c;
            % The above line does not cause the slider callback to be called.
            set_text_B2c()
            update()
        else
            set_text_B2c()
        end
    end
    function text_B2s_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_B2s.Max]),slider_B2s.Min]);
            B2s = val;
            slider_B2s.Value = B2s;
            % The above line does not cause the slider callback to be called.
            set_text_B2s()
            update()
        else
            set_text_B2s()
        end
    end
    function text_r0_callback(source,callbackdata)
        [val,count] = sscanf(source.String,'%g');
        if count==1
            val = max([min([val,slider_r0.Max]),slider_r0.Min]);
            r0 = val;
            slider_r0.Value = r0;
            % The above line does not cause the slider callback to be called.
            set_text_r0()
            update()
        else
            set_text_r0()
        end
    end

    function set_text_nfp
        text_nfp.String = sprintf('%g',nfp);
    end
    function set_text_R0c1
        text_R0c1.String = sprintf('%g',R0c1);
    end
    function set_text_Z0s1
        text_Z0s1.String = sprintf('%g',Z0s1);
    end
    function set_text_eta_bar
        text_eta_bar.String = sprintf('%g',eta_bar);
    end
    function set_text_sigma
        text_sigma.String = sprintf('%g',sigma);
    end
    function set_text_I2
        text_I2.String = sprintf('%g',I2);
    end
    function set_text_minus_mu0_p2
        text_minus_mu0_p2.String = sprintf('%g',minus_mu0_p2);
    end
    function set_text_B2c
        text_B2c.String = sprintf('%g',B2c);
    end
    function set_text_B2s
        text_B2s.String = sprintf('%g',B2s);
    end
    function set_text_r0
        text_r0.String = sprintf('%g',r0);
    end

    function update()
        % Write quasisymmetry input namelist:
        input_filename = [quasisymmetry_temp_files_directory,'quasisymmetry_in.matlab'];
        fid = fopen(input_filename,'w');
        if fid<0
            error(['Unable to open file ',input_filename,' for writing.'])
        end
        fprintf(fid,'&quasisymmetry\n');
        fprintf(fid,' general_option="single"\n');
        fprintf(fid,' untwist=.false.\n');
        if finite_r_nonlinear_button.Value == finite_r_nonlinear_button.Max
            fprintf(fid,' finite_r_option="nonlinear"\n');
        else
            fprintf(fid,' finite_r_option="linear"\n');
        end
        %{
        if order_r_squared_button.Value == order_r_squared_button.Max
            fprintf(fid,' order_r_squared=.true.\n');
        else
            fprintf(fid,' order_r_squared=.false.\n');
        end
        %}
        order_r_option = get(order_r_option_button_group,'SelectedObject');
        if order_r_option == order_r_option_r1_button
            fprintf(fid,' order_r_option="r1"\n');
        elseif order_r_option == order_r_option_r2_button
            fprintf(fid,' order_r_option="r2"\n');
        elseif order_r_option == order_r_option_r3_flux_constraint_button
            fprintf(fid,' order_r_option="r3_flux_constraint"\n');
        elseif order_r_option == order_r_option_r3_simplified_button
            fprintf(fid,' order_r_option="r3_simplified"\n');
        elseif order_r_option == order_r_option_r3_simplified_with_Z3_button
            fprintf(fid,' order_r_option="r3_simplified_with_Z3"\n');
        else
            error('Should not get here!')
        end
        fprintf(fid,' N_phi=101\n');
        fprintf(fid,' nfp = %d\n',round(nfp));
        fprintf(fid,' eta_bar = %.15g\n',eta_bar);
        fprintf(fid,' sigma_initial = %.15g\n',sigma);
        fprintf(fid,' I2_over_B0 = %.15g\n',I2);
        fprintf(fid,' R0c = 1, %.15g\n',R0c1);
        fprintf(fid,' Z0s = 0, %.15g\n',Z0s1);
        fprintf(fid,' B2s = %.15g\n',B2s);
        fprintf(fid,' B2c = %.15g\n',B2c);
        fprintf(fid,' r = %.15g\n',r0);
        fprintf(fid,[' vmec_template_filename = "',quasisymmetry_executable_directory,'input.li383_vacuum"\n']);
        fprintf(fid,'/\n');
        fclose(fid);
        
        % Run the executable:
        system(['cd ',quasisymmetry_temp_files_directory,';',quasisymmetry_executable_directory,'quasisymmetry quasisymmetry_in.matlab > output']);
        number_of_calls = number_of_calls + 1;
        %fprintf('%d\n',number_of_calls);
        
        output_filename = [quasisymmetry_temp_files_directory,'quasisymmetry_out.matlab.nc'];
        
        iota = ncread(output_filename,'iota');
        elongation = ncread(output_filename,'max_elongation');
        helicity = ncread(output_filename,'axis_helicity');
        %B20_mean = ncread(output_filename,'B20_mean');
        RBC = ncread(output_filename,'RBC');
        RBS = ncread(output_filename,'RBS');
        ZBC = ncread(output_filename,'ZBC');
        ZBS = ncread(output_filename,'ZBS');
        
        label_iota.String = sprintf('Iota: %.4g',iota);
        label_elongation.String = sprintf('Elongation: %.4g',elongation);
        if helicity==0
            label_helicity.String = sprintf('Quasi-axisymmetry');
        else
            label_helicity.String = sprintf('Quasi-helical symmetry');
        end
        %label_B20.String = sprintf('B20: %g',B20);
        
        R = zeros(size(phi2D));
        Z = zeros(size(phi2D));
        ntor = (size(RBC,1)-1)/2;
        mpol = size(RBC,2)-1;
        for m = 0:mpol
            for n = (-ntor):ntor
                angle = m * theta2D - nfp * n * phi2D;
                sinangle = sin(angle);
                cosangle = cos(angle);
                R = R + RBC(n+ntor+1,m+1) * cosangle + RBS(n+ntor+1,m+1) * sinangle;
                Z = Z + ZBC(n+ntor+1,m+1) * cosangle + ZBS(n+ntor+1,m+1) * sinangle;
            end
        end
        X = R .* cosphi;
        Y = R .* sinphi;
        
        B = 1 + eta_bar * r0 * cos(theta2D);
        
        phi_decimated = linspace(0,2*pi/nfp,5);
        phi_decimated(end) = [];
        [phi2D_decimated, theta2D_decimated] = meshgrid(phi_decimated,theta);
        R_decimated = zeros(size(phi2D_decimated));
        Z_decimated = zeros(size(phi2D_decimated));
        for m = 0:mpol
            for n = (-ntor):ntor
                angle = m * theta2D_decimated - nfp * n * phi2D_decimated;
                sinangle = sin(angle);
                cosangle = cos(angle);
                R_decimated = R_decimated + RBC(n+ntor+1,m+1) * cosangle + RBS(n+ntor+1,m+1) * sinangle;
                Z_decimated = Z_decimated + ZBC(n+ntor+1,m+1) * cosangle + ZBS(n+ntor+1,m+1) * sinangle;
            end
        end

        %R_decimated = R(:,1:Nphi_stride:(Nphi_stride*3+1));
        %Z_decimated = Z(:,1:Nphi_stride:(Nphi_stride*3+1));
        assignin('base','R_decimated',R_decimated)
        assignin('base','Z_decimated',Z_decimated)
        
        if first_update
            first_update = false;
            
            surf_handle = surf(X,Y,Z,B);
            daspect([1,1,1])
            axis vis3d off
            set(gca,'clipping','off','color','w')
            zoom(1.7)
            rotate3d on
            
            ax2 = axes('Units','pixels','Position',[260,25,200,200]);
            slices_plot_handle = plot(R_decimated,Z_decimated,'XDataSource','R_decimated','YDataSource','Z_decimated');
            xlim([0.8,1.2])
            axis equal
            %xlabel('R')
            %ylabel('Z')
        else
            set(surf_handle,'XData',X,'YData',Y,'ZData',Z,'CData',B);
            %set(slices_plot_handle,'XData',R_decimated,'YData',Z_decimated);

            refreshdata
        end
        
        
    end


end