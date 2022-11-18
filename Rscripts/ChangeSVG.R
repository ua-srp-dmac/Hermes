library(stringr)

ChangeSVGBasic <- function(filename1,temp){
  
  x<- readLines(filename1)
  x2 <- x;
  flag1 = 0;
  
  
  for(i in 1:length(x)){
    
    if(str_detect(x[i], "<svg")){
      
      str_sub(x2[i],1,4) <- '<svg id="tooltip-svg"'
      
    } else if(str_detect(x[i], "<text") && flag1 != 3){
      flag1 = 1
      res <- str_match(x[i], ">(.*?)<")
      x2[i] <- gsub('textLength=(.*?) ',paste0('class = "test1" data-tooltip-text="', temp[temp[,2] ==  res[2],1] , '" '),x[i])
      x2[i] <- gsub('_(.*?)<', "<",x2[i])
    } else if(str_detect(x[i], "</svg>")){
      
      x2[i+1] = "</svg>"
      x2[i] = 
        '
  <g id="tooltip" visibility="hidden">
    <rect x="2" y="2" width="80" height="24" style="fill: black; opacity:0.4;" rx="2" ry="2"/>
      <rect width="80" height="24" style="fill: white;" rx="2" ry="2"/>
        <text x="4" y="16">Tooltip</text>
          </g>
          
          
          <script type="text/ecmascript"><![CDATA[
            (function() {
              var svg = document.getElementById("tooltip-svg");
              var tooltip = svg.getElementById("tooltip");
              var tooltipText = tooltip.getElementsByTagName("text")[0];
              var tooltipRects = tooltip.getElementsByTagName("rect");
              var triggers = svg.getElementsByClassName("test1");
              for (var i = 0; i < triggers.length; i++) {
                triggers[i].addEventListener("mousemove", showTooltip);
                triggers[i].addEventListener("mouseout", hideTooltip);
              }
              function showTooltip(evt) {
                var CTM = svg.getScreenCTM();
                var x = (evt.clientX - CTM.e + 6) / CTM.a;
                var y = (evt.clientY - CTM.f + 10) / CTM.d;
                tooltip.setAttributeNS(null, "transform", "translate(" + x + " " + y + ")");
                tooltip.setAttributeNS(null, "visibility", "visible");
                tooltipText.firstChild.data = evt.target.getAttributeNS(null, "data-tooltip-text");
                var length = tooltipText.getComputedTextLength();
                for (var i = 0; i < tooltipRects.length; i++) {
                  tooltipRects[i].setAttributeNS(null, "width", length + 8);
                }
              }
              function hideTooltip(evt) {
                tooltip.setAttributeNS(null, "visibility", "hidden");
              }
            })()
            ]]></script>
  '
      
    } else{
      
      if(flag1 == 1){
        flag1 = 3
      }
    }
    
  }
  
  write(x2, filename1)
  
}

ChangeSVGUPSET <- function(filename1,temp){
  
  x<- readLines(filename1)
  x2 <- x;
  flag1 = 0;
  
  
  for(i in 1:length(x)){
    
    if(str_detect(x[i], "<svg")){
      
      str_sub(x2[i],1,4) <- '<svg id="tooltip-svg"'
      
    } else if(str_detect(x[i], "<text")){
      if(flag1 == 0){
	flag1 = 1
	} else if(flag1 == 3){
	
      res <- str_match(x[i], ">(.*?)<")
      x2[i] <- gsub('textLength=(.*?) ',paste0('class = "test1" data-tooltip-text="', temp[temp[,2] ==  res[2],1] , '" '),x[i])
      x2[i] <- gsub('_(.*?)<', "<",x2[i])
	}
    } else if(str_detect(x[i], "</svg>")){
      
      x2[i+1] = "</svg>"
      x2[i] = 
        '
  <g id="tooltip" visibility="hidden">
    <rect x="2" y="2" width="80" height="24" style="fill: black; opacity:0.4;" rx="2" ry="2"/>
      <rect width="80" height="24" style="fill: white;" rx="2" ry="2"/>
        <text x="4" y="16">Tooltip</text>
          </g>
          
          
          <script type="text/ecmascript"><![CDATA[
            (function() {
              var svg = document.getElementById("tooltip-svg");
              var tooltip = svg.getElementById("tooltip");
              var tooltipText = tooltip.getElementsByTagName("text")[0];
              var tooltipRects = tooltip.getElementsByTagName("rect");
              var triggers = svg.getElementsByClassName("test1");
              for (var i = 0; i < triggers.length; i++) {
                triggers[i].addEventListener("mousemove", showTooltip);
                triggers[i].addEventListener("mouseout", hideTooltip);
              }
              function showTooltip(evt) {
                var CTM = svg.getScreenCTM();
                var x = (evt.clientX - CTM.e + 6) / CTM.a;
                var y = (evt.clientY - CTM.f + 10) / CTM.d;
                tooltip.setAttributeNS(null, "transform", "translate(" + x + " " + y + ")");
                tooltip.setAttributeNS(null, "visibility", "visible");
                tooltipText.firstChild.data = evt.target.getAttributeNS(null, "data-tooltip-text");
                var length = tooltipText.getComputedTextLength();
                for (var i = 0; i < tooltipRects.length; i++) {
                  tooltipRects[i].setAttributeNS(null, "width", length + 8);
                }
              }
              function hideTooltip(evt) {
                tooltip.setAttributeNS(null, "visibility", "hidden");
              }
            })()
            ]]></script>
  '
      
    } else{
      
      if(flag1 == 1){
        flag1 = 3
      }
    }
    
  }
  
  write(x2, filename1)
  
}
