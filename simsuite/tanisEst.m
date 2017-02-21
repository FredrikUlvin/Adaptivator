classdef tanisEst < handle
    properties
        innoWindow = [];
        covarWindow = 0;
        
        sigLvl = 0.05;
        chi2Lim = [];
        tanis = 0;
        
        alpha = 1;
        tuneStep = 2;
        counter = 0;
    end
    
    methods 
        function self = tanisEst(measSize, alphaStep, sampleSize)
            self.innoWindow = zeros(measSize, sampleSize);
            self.updateChiLim(measSize, sampleSize);
            self.counter = sampleSize;
            self.tuneStep = alphaStep;
        end
        
        function alpha = update(self, inno, S)
            self.innoWindow = [inno self.innoWindow(:,1:end-1)];
            
            if self.counter > 0
                self.counter = self.counter - 1;
            else
                self.tanisCalc(S);
                self.tuneAlpha(self.tanis);
                self.covarWindowCalc();
            end
            alpha = self.alpha;
        end
        
        function covarWindowCalc(self)           
            tempc = 0;
            for i = 1:length(self.innoWindow)
                tempc = tempc + self.innoWindow(:,i)' * self.innoWindow(:,i);
            end
            self.covarWindow = (1/length(self.innoWindow)) * tempc;
        end
        
        function nis = nisCalc(self, inno, S)
            nis = inno' * inv(S) * inno;
        end
        
        function tanisCalc(self, S)
            temp = 0;
            for i = 1:length(self.innoWindow)
                temp = temp + self.innoWindow(:,i)' * inv(S) * self.innoWindow(:,i);
            end
            self.tanis = temp/length(self.innoWindow(1,:));
        end
        
        function updateChiLim(self, meas, sample)
            self.chi2Lim =  (1/sample) .* chi2inv([self.sigLvl (1 - self.sigLvl)] , sample*meas);
        end
        
        function tuneAlpha(self, test)
            if test < self.chi2Lim(1)
                if self.alpha > self.tuneStep^(-10)
                    self.alpha = self.alpha / self.tuneStep;
                else
                    self.alpha = self.alpha;
                end
            elseif test > self.chi2Lim(2)
                if self.alpha < self.tuneStep^10
                    self.alpha = self.alpha * self.tuneStep;
                else
                    self.alpha = self.alpha;
                end
            else
                self.alpha = self.alpha;
            end
        end
    end
end