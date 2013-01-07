function func = getpoints
    f = figure;
    axis image;
    axis([0 100 0 100]);
    axis manual;
    hold on;
    a = gca;
    points = zeros(0, 3);
    function innerPoints = getPointsFunc()
        innerPoints = points;
    end
    func = @getPointsFunc;
    function clickFunc(~, ~, ~)
        fprintf('clickFunc called\n');
        button = get(f, 'SelectionType');
        point = get(a, 'CurrentPoint');
        point = point(1,1:2);
        cls = -1;
        ch = 'xr';
        if ~isequal(button, 'normal')
            cls = 1;
            ch = 'ob';
        end
        points = [points; cls, point];
        plot(point(1), point(2), ch);
    end
    set(a, 'ButtonDownFcn', @clickFunc);
end