import QtQuick 2.4
import QtQuick.Window 2.2
import "Creat.js" as Monitor

Window {
    visible: true
    width: 400
    height: 400
    Rectangle{
        id:rec
        x:50
        y:50
        width:250
        height: 250
        color: "gray"

    Canvas{
       x:50
       y:50
       width: 150
       height: 150
    contextType: "2d"
    visible: true
    onPaint: {//绘图事件的响应
        context.lineWidth=20;
        context.strokeStyle="red";
        context.fillStyle="blue";
        context.beginPath();
        context.rect(10,0,120,120);//一定要注意画图的区域不能超过画布的大小，不然会看不到或者只看到一部分
        context.fill();
        context.stroke();


    }


    }
    }
}
