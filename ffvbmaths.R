## FFVB definitions for this model



$$\pi(y_{i} | p, \tau) = \frac{1}{\sqrt{(2\pi)}} (\sum_{k=1}^Kp_k^2\sigma_{k}^2 + \tau^{-1})^{\frac{n}{2}} \exp\left(-\frac{(\sum_{k=1}^Kp_k^2\sigma_{k}^2 + \tau^{-1})}{2} \sum_{i=1}^n\left(y_{i} - \sum_{k=1}^Kp_k\mu_{k}\right)^2\right)$$
  
  
  
  
  $$h(\theta) = \log(\pi(p)\pi(\tau)\pi(y|p, \tau))$$
    
    
    $\log{q_\lambda}(\theta) = \log(q(p)q(\tau))$
      
      
      $$h(\theta) = -\log(\beta(\alpha_0)) + \sum_{k=1}^K\left((\alpha_0-1)(\log(p_k))\right) + c_0\log(d_0) - \log(\Gamma(c_0)) +
        (c_0-1)\log(\tau) - d_0\tau$$
        $$-\frac{1}{2}\log(2\pi) +\frac{n}{2}\log(\tau) - \frac{\tau}{2}\left(\sum_{i=1}^n\left(y_{i} - \sum_{k=1}^K(p_k\mu_{k})\right)^2\right)$$
        
        
        
        
        $$\log(q_\lambda(\theta)) = -\log(\beta(\alpha)) + \sum_{k=1}^K\left((\alpha-1)(\log(p_k))\right)$$ 
          
          
          Finally, the derivative of the variational approximation for each of the four parameters is given by:
          
          wrt c
        $$\log(d) - \frac{\Gamma'(c)}{\Gamma(c)} + log(\tau_j)$$

wrt d
$$\frac{c}{d} - \tau_j$$

wrt each alpha
$$\log(p_1) - \frac{1}{\beta(\alpha)}*\beta(\alpha)\left(\Psi(\alpha_1)-\Psi(\sum_{k=1}^K(\alpha_k)\right)$$