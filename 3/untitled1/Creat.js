.pragma library

function clear(){
}
//划线
function drawLine(pstart, pend, opts){
    var ph = 1;
    var pw = 1;
    var color = "DarkRed";
    if(opts){
        color = opts.color ? opts.color: color;
    }
    var slope;
    var noSlope = false;
    var hdist = pend[0] - pstart[0];
    var vdist = pend[1] - pstart[1];
    if(hdist != 0){
        slope =  Math.abs(vdist/hdist);
    }else{
        noSlope = true;
    }
    var gapp = pw > ph ? ph : pw;

    var diagonal = Math.sqrt(Math.pow(hdist,2) + Math.pow(vdist,2));
    var pn = parseInt(diagonal/gapp);
    if(pn < 3){pn=pn*3+1};
    var vgap = Math.abs(vdist)/pn;
    var hgap = Math.abs(hdist)/pn;
    for(var i = 0; i< pn ; i++){
   //画点
        drawPoint({
            pw: pw,
            ph: ph,
            color: color,
            point: [(hgap*i*(pend[0]<pstart[0]?-1:1)*(noSlope?0:1)+pstart[0]),(vgap*i*(pend[1]<pstart[1]?-1:1)+pstart[1])]
        });
    }
}

function drawPoint(opts){
    document.write("<span id='"+opts.point[0]+""+opts.point[1]+"' style='display: inline-block; width: "+opts.pw+"px; height: "+opts.ph+"px; background-color: "+opts.color+"; position: absolute "+opts.point[0]+"px; top: "+opts.point[1]+"px;'>"+(opts.point[2]?("<div style='position: relative;'><span style='position: absolute; left: 5px; bottom: 1px; text-align: left; width: 100px;'>"+opts.point[2]+"</span></div>"):"")+"</span>");
}
//画矩形
function drawRectangle(leftTop, width, high){
    drawPolygon([
        leftTop,
        [leftTop[0], leftTop[1]+high],
        [leftTop[0]+width, leftTop[1]+high],
        [leftTop[0]+width, leftTop[1]]
    ]);
    //填充
    //document.write("<span style='height: "+(high-1)+"px; width: "+(width-1)+"px; background-color: "+"Green"+"; position: absolute; left:"+(leftTop[0]+1)+"px; top: "+(leftTop[1]+1)+"'></span>");
}
