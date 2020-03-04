% This page calculates pq loads
classdef pqloads<matlab.mixin.Copyable    
    properties
        mpf;                % basic model 
        qlist;              % generated mpflist
        load;               % load csv
        bus;                % bus num
        
        halfhload;          %half hour load
        halftlist;          % 1:0.5:24
        interploadlist;     % interplot p
        interqloadlist;     % q
        interptlist;        %t
        
        p0list;             %mpf.bus 3
        q0list;             %mpf.bus 4
        
        qstlist;
        
        isrenew;
        renewables;         % new energy only
        loc;                % renewable location(bus)
        inttlist;
        err;                % additional error
        tol;                % tol for ac power flow calculation
        
    end
    properties(Dependent)
        Modulus
    end
    
    methods
        function obj = pqloads(mpf,varargin)% load dt renewables loc
            obj.mpf = mpf;      % load model mpf 
            obj.pindex;         % calculate p0 q0
            obj.err = 0.001;    %default err
            obj.tol = 10^-8;
            obj.loc;
            obj.inttlist = 0:1:24;
            
            obj.isrenew = 0;
            obj.load = varargin{1};   
            if nargin == 6 
            obj.renewables = varargin{3};
            obj.loc = varargin{4};
            obj.err = varargin{5};
            obj.isrenew = 1;
            elseif nargin == 5
            obj.err = varargin{3}; 
            obj.tol = varargin{4};
            end
            obj.halftlist = 0:0.5:24; 
            obj.halfhourload;   % build standard load
            
            obj.interpload(varargin{2});
            
            obj.buildqlist;
            
        end
        
        function obj = halfhourload(obj)
            [m,~] = size(obj.load);
            cutlist = fix(linspace(1,m,obj.bus+1));
            loadltp = obj.load;
            loadbus = [];
            for i  = 1:obj.bus
                s = cutlist(i);
                e = cutlist(i+1);
                loadtp = loadltp(s:e,:);
                loadtpsum = sum(loadtp);
                loadbus = [loadbus; loadtpsum/max(loadtpsum)];
            end
            obj.halfhload = loadbus;             
        end
        
        function plothload(obj)
            figure
            plot(transpose(obj.halfhload))
        end
        
        function obj = interpload(obj,dt)
            interplist = 0:dt:24-dt;
            obj.interptlist = interplist;
            interloadlist = [];
            interloadqlist = [];
            htp = obj.halfhload;
            lna = obj.renewables;
            loclist = obj.loc;
            u = 1;
            for i  = 1:obj.bus              
                intp0 = interp1(obj.halftlist,htp(i,:),interplist);
                intp = intp0*0.3;
                intq = intp0*0.12;
             
                [li,~] = ismember(i, loclist);
                intprenew = [];
                if li == 1 %exist
                    intprenew = interp1(0:24,lna(u,:),interplist)*0.2;
                    u = u + 1;
                end
                
                for k = 1:length(intp)
                    intp(k) = intp(k)*(1+ 0.1*randn());
                    if li == 1
                    intp(k) = intp(k)- intprenew(k);
                    end
                end
                
                for k = 1:length(intq)
                    intq(k) = intq(k)*(1+ 0.1*randn());
                    if li == 1
                    intq(k) = intq(k)- intprenew(k)*0.5;
                    end
                end
                interloadlist = [interloadlist; intp];
                interloadqlist = [interloadqlist; intq];
            end
            obj.interploadlist = interloadlist;
            obj.interqloadlist = interloadqlist;
        end
        
      
        
        function plotinterpload(obj,varargin)
            figure
            hold on;
            if nargin == 1
                 plot(transpose(obj.interploadlist));
            else 
                 plot(obj.interptlist,transpose(obj.interploadlist(varargin{1},:)));
                 plot(obj.interptlist,transpose(obj.interqloadlist(varargin{1},:)));
                 grid on
                 xlim([0 24]);
                 set(gca,'XTick',0:3:24);
                 xlabel('Time/Hour');
                 ylabel({'Power/MW (MVA)'});
                 title(strcat('Bus:',num2str(varargin{1}+1)));
                 legend('p','q')
                 ylim([ 0 0.4]);
                
            end
%             saveas(gcf,['renewload\',num2str(varargin{1}+1),'.eps'],'psc2')
%             saveas(gcf,['renewload\',num2str(varargin{1}+1),'.png'])
        end
        
        function saveloadeps(obj)
            for i = 1:obj.bus
                obj.plotinterpload(i);
            end
        end
        
        function obj = pindex(obj)
            bus = obj.mpf.bus;
            obj.p0list = bus(:,3);
            obj.q0list = bus(:,4);
            obj.bus = length(obj.p0list)-1;
        end
        function b = erb(obj,a,thro)
            b = a;
              if a>thro
                  b = thro;
              end
              if a<-thro
                  b = -thro;
              end
        end
        
        function obj = buildqlist(obj)
            
            qlist = [];
            qstlist = [];
            for i = 1:length(obj.interptlist)
                if mod(i,50) == 0 ||i == length(obj.interptlist)
                    fprintf('[data processing] Build No. %d dataset\n',i)
                end
                mpf1 = obj.mpf;
                mpf1.bus(2:obj.bus+1,3) = obj.interploadlist(:,i);
                mpf1.bus(2:obj.bus+1,4) = obj.interqloadlist(:,i);
                mpopt = mpoption('verbose',0,'out.all',0);
                q_stand = (runpf(mpf1,mpopt)); %ฑ๊ปฏ
                
                
                mpopt.pf.tol = obj.tol;
                
                q = (runpf(mpf1,mpopt));
                
                for k = 2:obj.bus+1
                    q.bus(k,3) =  q.bus(k,3)*(1+obj.err*randn());
                     
                end
                
                for k = 2:obj.bus+1
                   q.bus(k,4) =  q.bus(k,4)*(1+obj.err*randn());
                end
%                 for k = 2:obj.bus+1
%                    q.bus(k,8) =  q.bus(k,8)+0.00001*randn();
%                 end
% %                 
                qlist = [qlist q];
                qstlist = [qstlist q_stand];
            end
            obj.qlist = qlist;
            obj.qstlist = qstlist;
        end
        
    end

end
