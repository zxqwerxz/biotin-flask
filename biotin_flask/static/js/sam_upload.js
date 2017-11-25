$(function() {
  $('#samfiles tbody button').click(function(event) {
    $filename = $(this).attr('title');
    data = { 'file' : $filename };
    var that = $(this)
    $.ajax({
      type: 'POST',
      url: $SCRIPT_ROOT + '/sam/delete',
      data: JSON.stringify(data),
      contentType: 'application/json',
      success: function(result) {
        that.closest('tr').remove();
        if ($('#samfiles tbody tr').length == 0) {
          $('#samfiles tbody').append('<tr><td>No files currently uploaded.</td></tr>');
          $('#samfiles #deleteAll').remove();
        }
      }
    });
  });

  $('#samfiles #deleteAll').click(function(event) {
    $.ajax({
      type: 'POST',
      url: $SCRIPT_ROOT + '/sam/delete_allbam',
      data: JSON.stringify({ }),
      contentType: 'application/json',
      success: function(result) {
        $('#samfiles tbody tr').remove();
        if ($('#samfiles tbody tr').length == 0) {
          $('#samfiles tbody').append('<tr><td>No files currently uploaded.</td></tr>');
          $('#samfiles #deleteAll').remove();
        }
      }
    });
  });

  $('#fastafile tbody button').click(function(event) {
    $filename = $(this).attr('title');
    data = { 'file' : $filename };
    var that = $(this)
    $.ajax({
      type: 'POST',
      url: $SCRIPT_ROOT + '/sam/delete',
      data: JSON.stringify(data),
      contentType: 'application/json',
      success: function(result) {
        that.closest('tr').remove();
        if ($('#fastafile tbody tr').length == 0) {
          $('#fastafile tbody').append('<tr><td>No file currently uploaded.</td></tr>');
          $('#fastaupload').css('display', 'block');
        }
      }
    });
  });
});
