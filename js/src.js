t_max = 1000;
t = 0;
tAutoUp = true;
mouseDown = false;
startX = null;
startT = null;

var dom = document.getElementById('chart-container');
var myChart = echarts.init(dom, null, {
  renderer: 'canvas',
  useDirtyRect: false
});
var app = {};

var option;

option = {
  title: {
    text: 't = 0'
  },
  tooltip: {
    trigger: 'axis',
    axisPointer: {
      type: 'cross',
      label: {
        backgroundColor: '#6a7985'
      }
    }
  },
  legend: {
    data: ['模态叠加法', 'Newmark-beta', 'Wilson-theta',]
  },
  grid: {
    left: '3%',
    right: '4%',
    bottom: '3%',
    containLabel: true
  },
  xAxis: [
    {
      min: 0.0,
      max: 1.0,
    }
  ],
  yAxis: [
    {
      min: -0.1,
      max: 0.1,
    }
  ],
  series: [
    {
      name: '模态叠加法',
      type: 'line',
      symbol: "none",
      animation: false,
      areaStyle: {},
      emphasis: {
        focus: 'series'
      },
      data: data1[t]
    },
    {
      name: 'Newmark-beta',
      type: 'line',
      symbol: "none",
      animation: false,
      areaStyle: {},
      emphasis: {
        focus: 'series'
      },
      data: data2[t]
    },
    {
      name: 'Wilson-theta',
      type: 'line',
      symbol: "none",
      animation: false,
      areaStyle: {},
      emphasis: {
        focus: 'series'
      },
      data: data3[t]
    },
  ]
};

myChart.getZr().on('click', function(params) {
  if(tAutoUp){
    tAutoUp = false;
  }
  else{
    tAutoUp = true;
  }
});

myChart.getZr().on('mousedown', function(params) {
  mouseDown = true;
  startX = params.offsetX;
  startT = t;
});

myChart.getZr().on('mouseup', function(params) {
  mouseDown = false;
  startX = null;
  startT = null;
});

myChart.getZr().on('mousemove', function(params) {
  if(mouseDown){
    t = startT + parseInt((params.offsetX - startX));
  }
  if(t < 0){
    t = 0;
  }
  if(t > t_max){
    t = t_max;
  }
});

var int=self.setInterval("clock()", 20);
function clock(){
  if(tAutoUp && !mouseDown){
    if(t < t_max){
      t += 1;
    }
  }
  myChart.setOption(
    {title: {
      text: 't = ' + t/(10*t_max),
      },
      
      series: [
      {
        name: '模态叠加法',
        type: 'line',
        symbol: "none",
        animation: false,
        areaStyle: {},
        emphasis: {
          focus: 'series'
        },
        data: data1[t]
      },
      {
        name: 'Newmark-beta',
        type: 'line',
        symbol: "none",
        animation: false,
        areaStyle: {},
        emphasis: {
          focus: 'series'
        },
        data: data2[t]
      },
      {
        name: 'Wilson-theta',
        type: 'line',
        symbol: "none",
        animation: false,
        areaStyle: {},
        emphasis: {
          focus: 'series'
        },
        data: data3[t]
      },
    ]
  }
  )
} 

myChart.setOption(option);

window.addEventListener('resize', myChart.resize);