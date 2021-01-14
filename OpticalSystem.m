classdef OpticalSystem < handle
    properties
        Elements = {}; % Optical Elements
        Locs = []; % Location of optical elements
        StartPoint = 0; % location of starting point
        unit = 'mm'; % string for unit
        ElemDiam = []; % Optical Element Diameter for plotting purposes
    end
    methods
        function self = OpticalSystem(sp, unit)
            self.StartPoint = sp;
            if exist('unit','var')
                self.unit = unit;
            end
        end
        function self = setStartPoint(self, sp)
            self.StartPoint = sp;
        end
        function self = addLens(self, z, f, diam)
            if z < self.StartPoint
                error('Place the Lens after the StartPoint');
            end
            self.Locs(end + 1) = z;
            self.Elements{end + 1} = self.Lens(f);
            if ~exist('diam','var')
                diam = 25;
            end
            self.ElemDiam(end + 1) = diam;
            self.sortLocs();
        end
        function res = propagateRays(self, rays, zlist, colors)
            % rays should be a 2 x N matrix, Zlist is a list of points to
            % propagate to and plot to 
            temp = size(rays);
            N = temp(2); % number of columns
            res = zeros(2, N, length(zlist));
            for i = 1:N
                for j = 1:length(zlist)
                    thisray = rays(:,i);
                    thisr = zlist(j);
                    rstart = self.StartPoint;
                    if thisr < rstart
                        error('Choose r that is beyond starting point')
                    end
                    passed_elements = find(thisr > self.Locs);
                    for k = passed_elements
                        thisray = self.FreeProp(self.Locs(k) - rstart) * thisray; % free propagate to element
                        thisray = self.Elements{k} * thisray; %apply element
                        rstart = self.Locs(k); % advance rstart
                    end
                    % final free prop
                    thisray = self.FreeProp(thisr - rstart) * thisray;
                    res(:,i,j) = thisray;
                end
            end
            %plotting
            figure();
            hold on;
            for i = 1:N
                if exist('colors','var')
                    plot(zlist, reshape(res(1,i,:),[1,length(zlist)]),colors{i})
                else
                    plot(zlist, reshape(res(1,i,:),[1,length(zlist)]))
                end
            end
            plot(zlist, zeros(1,length(zlist)), 'k--');
            for i = 1:length(self.ElemDiam)
                plot([1,1] * self.Locs(i), self.ElemDiam(i) * [1/2, -1/2], 'k-');
            end
            xlabel(sprintf('z Position (%s)', self.unit))
            ylabel(sprintf('Distance From Axis (%s)', self.unit))
        end
    end
    methods (Static)
        function res = Lens(f)
            res = [[1,0];[-1/f,1]];
        end
        function res = FreeProp(L)
            res = [[1, L];[0,1]];
        end
    end
    methods(Access = private)
        function self = sortLocs(self)
            [sorted, ind] = sort(self.Locs);
            self.Locs = sorted;
            self.Elements = self.Elements(ind);
            self.ElemDiam = self.ElemDiam(ind);
        end
    end
    
end